import psycopg2
import argparse
import logging
import time

from google.protobuf.json_format import Parse, MessageToJson
from google.protobuf.timestamp_pb2 import Timestamp

from phenopackets import Phenopacket,Individual, Sex, PhenotypicFeature, OntologyClass, \
		TimeElement,VitalStatus,Quantity,Measurement,Value

from datetime import datetime


def get_postgres_connection(dbname,username,password,host,port):
	return psycopg2.connect(dbname=dbname,user=username,password=password,host=host,port=port)


# sql queries

def get_individual_query(patient_id):
	return f"SELECT p.person_id::text as id, \
		null as alternate_ids,\
		p.birth_datetime as date_of_birth,\
		max(vo.visit_start_date) as time_at_last_encounter,\
		(case when d.person_id is null then 0 else 2 end) as vital_status, \
		(case when p.gender_concept_id is null then 0\
		 when p.gender_concept_id = 8532 then 1\
		 when p.gender_concept_id = 8507 then 2\
		 else 3 end) as sex,\
		null as karyotypic_sex,\
		null as gender,\
		'NCBITaxon:9606' as taxonomy_id,\
		'human' as taxonomy_label\
		FROM cdm_synthea10.person p\
		LEFT JOIN cdm_synthea10.visit_occurrence vo on vo.person_id = p.person_id\
		LEFT JOIN cdm_synthea10.death d on d.person_id = p.person_id\
		WHERE p.person_id={patient_id} \
		GROUP BY p.person_id, p.birth_datetime, vital_status, sex\
		ORDER BY p.person_id asc ;" \

def get_vitalstatus_query(patient_id):
	return f"SELECT d.person_id, \
		(case when d.person_id is null then 0 else 2 end) as vital_status, \
		d.death_datetime as time_of_death, \
		(case when c.concept_code is null then null else concat(c.vocabulary_id,':',c.concept_code) end) as cause_of_death_id, \
		c.concept_name as cause_of_death_label\
		FROM cdm_synthea10.person p\
		LEFT JOIN cdm_synthea10.death d on d.person_id = p.person_id \
		LEFT JOIN cdm_synthea10.condition_occurrence co on co.person_id = p.person_id and co.condition_concept_id = d.cause_concept_id \
		LEFT JOIN cdm_synthea10.concept c on c.concept_id = d.cause_concept_id\
		WHERE d.person_id={patient_id} \
		GROUP BY d.person_id, vital_status, time_of_death, cause_of_death_id, cause_of_death_label \
		ORDER BY d.person_id asc;" \

def get_measurement_query(patient_id):
	return f"select m.person_id,\
		m.measurement_concept_id,\
		concat(c.vocabulary_id,':',c.concept_code) as assay_id,\
		c.concept_name as assay_label,\
		m.value_as_number,\
		concat(c3.vocabulary_id,':',c3.concept_code) as value_id,\
		m.value_source_value as value_label,\
		m.range_low,\
		m.range_high,\
		m.measurement_datetime,\
		concat(c2.vocabulary_id,':',c2.concept_code) as unit_id, \
		c2.concept_name as unit_label,\
		c2.concept_id,\
		m.unit_source_value,\
		m.visit_occurrence_id,\
		row_number() over (partition by m.person_id, m.measurement_datetime, m.visit_occurrence_id)\
		FROM cdm_synthea10.measurement m \
		left join cdm_synthea10.concept c on c.concept_id = m.measurement_concept_id\
		left join cdm_synthea10.concept c2 on c2.concept_id = m.unit_concept_id\
		left join cdm_synthea10.concept c3 on c3.concept_id = m.value_as_concept_id\
		where m.person_id={patient_id};"




#parsing

def parse_Individual(records):
	fields=["id","alternate_ids","date_of_birth","time_at_last_encounter","vital_status","sex","karyotypic_sex","gender","taxonomy_id","taxonomy_label"]
	return {i:j for i,j in zip(fields,records[0]) if j is not None }

def parse_VitalStatus(records):
	fields=["person_id","vital_status","time_of_death","cause_of_death_id","cause_of_death_label"]
	return {i:j for i,j in zip(fields,records[0]) if j is not None }

def parse_Measurements(records):
	fields=["person_id","measurement_concept_id","assay_id","assay_label","value_as_number","value_id","value_label","range_low","range_high","measurement_datetime","unit_id","unit_label","concept_id","unit_source_value","visit_occurrence_id","row_number"]
	for r in records:
		logging.debug(f'RECORD {r}')
	measurements=[]
	values_nono=[None]
#	values_nono=[None,"None:No matching concept","No matching concept"]
	for r in records:
		measurements.append({i:j for i,j in zip(fields,r) if j not in values_nono })

	return measurements

#creating dict

def createDictIndividual(mydict,vsdict):
	#GENDER,KARYOTIPIC SEX,SURVIVAL TIME IN DAYS NOT MAPPED SO MISSING
	idict={}
	if('id' in mydict):
		idict['id']=mydict['id']
	if('alternate_ids' in mydict):
		idict['alternate_ids']=mydict['alternate_ids']
	if('date_of_birth' in mydict):
		idict['date_of_birth']=convert_time(mydict['date_of_birth'])
	if('time_at_last_encounter' in mydict):
		idict['time_at_last_encounter']=convert_time(mydict['time_at_last_encounter'])
	if('sex' in mydict):
		if(mydict['sex']==0):
			idict['sex']='UNKNOWN_SEX'
		elif(mydict['sex']==1):
			idict['sex']='FEMALE'
		elif(mydict['sex']==2):
			idict['sex']='MALE'
		elif(mydict['sex']==3):
			idict['sex']='OTHER_SEX'
	if('taxonomy_id' in mydict):
		idict['taxonomy']={'id':mydict['taxonomy_id'],'label':mydict['taxonomy_label']} 

	if('vital_status' in mydict and mydict['vital_status']!=0):
		tempdict={}
		if('status' in vsdict):
			tempdict={}
			if(vsdict['status']==0):
				tempdict['status']='UNKNOWN_STATUS'
			elif(vsdict['status']==1):
				tempdict['status']='ALIVE'
			elif(vsdict['status']==2):
				tempdict['status']='DECEASED'            
		if('time_of_death' in vsdict):
			tempdict['time_of_death']=convert_time(vsdict['time_of_death'])
		if('cause_of_death_id' in vsdict):
			tempdict['cause_of_death']={'id':vsdict['cause_of_death_id'],'label':vsdict['cause_of_death_label']}

		idict['vital_status']=tempdict

	return idict


def createListDictMeasurements(md):
#	"person_id","measurement_concept_id","assay_id","assay_label","value_as_number","value_id","value_label","range_low","range_high","measurement_datetime","unit_id","unit_label","concept_id","unit_source_value","visit_occurrence_id","row_number"	
	ilist=[]
	for m in md:
		logging.debug(f'Measurement: {m}')

		#discard if value_as_number and value_label are missing
		if(not any(k in m for k in ('value_as_number','value_label'))):
			logging.debug(f'row discarded')
			ilist.append({'discarded':'yes'})
			continue

		tempdict={}
		#assay
		tempdict['assay']={'id':m['assay_id'],'label':m['assay_label']}

		if('measurement_datetime' in m):
			tempdict['time_observed']=convert_time(m['measurement_datetime'])

		if('value_as_number' in m):#Measurement with quantity
			tdict={}
#			tdict['value']=m['value_as_number']
			if('.' in m['value_label']):
				tdict['value']=float(m['value_label'])
			else:
				tdict['value']=int(m['value_label'])
			tdict['unit']={'id':m['unit_id'],'label':m['unit_label']}
			tempdict['value']={'quantity':tdict}


			if(any(k in m for k in ('range_low','range_high'))):
				tdict={}
				tdict['unit']={'id':m['unit_id'],'label':m['unit_label']}
				if('range_low' in m):
					tdict['low']=m['range_low']
				if('range_high' in m):
					tdict['high']=m['range_high']
				tempdict['reference_range']=tdict				
		else:#Measurement with ontology
			tempdict['value']={'id':m['value_id'],'label':m['value_label']}

		ilist.append(tempdict)
	return ilist


#creating pheno
  
def createPhenoIndividual(individualdict):
#individual part	
	if('date_of_birth' in individualdict):
		individualdict['date_of_birth']=Timestamp(seconds=convert_time_toseconds(individualdict['date_of_birth']))
	if('time_at_last_encounter' in individualdict):
		individualdict['time_at_last_encounter']=TimeElement(timestamp=Timestamp(seconds=convert_time_toseconds(individualdict['time_at_last_encounter'])))
	if('taxonomy' in individualdict):
		tx=OntologyClass(id=individualdict['taxonomy']['id'],label=individualdict['taxonomy']['label'])
		individualdict['taxonomy']=tx	
#vital status part
	if('vital_status' in individualdict):
		if('time_of_death' in individualdict['vital_status']):
			individualdict['vital_status']['time_of_death']=TimeElement(timestamp=Timestamp(seconds=convert_time_toseconds(individualdict['vital_status']['time_of_death'])))
		if('cause_of_death' in individualdict['vital_status']):
			cd=OntologyClass(id=individualdict['vital_status']['cause_of_death']['id'], label=individualdict['vital_status']['cause_of_death']['label'])
			individualdict['vital_status']['cause_of_death']=cd
		vs=VitalStatus(**individualdict['vital_status'])
		individualdict['vital_status']=vs
	
	logging.debug(f'{individualdict}')
	return Individual(**individualdict)


def createPhenoMeasurement(ilist):
	measurements=[]
	for i in ilist:
		if('discarded' in i):
			continue
		#assay
		i['assay']=OntologyClass(id=i['assay']['id'],label=i['assay']['label'])
		if('id' in i['value']):#ontology
			i['value']=Value(ontology_class=OntologyClass(id=i['value']['id'],label=i['value']['label']))
		else:
			i['value']=Value(quantity=Quantity(unit=OntologyClass(id=i['value']['quantity']['unit']['id'], \
				                                   label=i['value']['quantity']['unit']['label']), \
			                    value=i['value']['quantity']['value']))
		if('time_observed' in i):
			i['time_observed']=TimeElement(timestamp=Timestamp(seconds=convert_time_toseconds(i['time_observed'])))

		measurements.append(Measurement(**i))

	return measurements


def createMetadata(myname):
	metadata={}
	metadata['created']=Timestamp(seconds=int(time.time()))
	metadata['created_by']=myname
	mdr=[]
	md={}
	md['id']='hp'
	md['name']=	'human phenotype ontology'
	md['url']='http://purl.obolibrary.org/obo/hp.owl'
	md['version']='2018-03-08'
	md['namespace_prefix']='HP'
	md['iri_prefix']='hp'
	mdr.append(md)
	metadata['resources']= mdr     
	metadata['phenopacket_schema_version']='2.0'
	logging.debug(f'metadata: {metadata}')
	return metadata

def createPheno(myid,meta_data,subject=None,phenotypic_features=None,measurements=None,biosamples=None,interpretations=None,diseases=None,medical_actions=None,files=None):
	pheno={}
	pheno['id']=myid
	if(subject != None):
		pheno['subject']=subject
	if(phenotypic_features!= None):
		pheno['phenotypic_features']=phenotypic_features
	if(measurements != None):
		pheno['measurements']=measurements
	if(biosamples != None):
		pheno['biosamples']=biosamples
	if(interpretations != None):
		pheno['interpretations']=interpretations
	if(diseases != None):
		pheno['diseases']=diseases
	if(medical_actions != None):
		pheno['medical_actions']=medical_actions
	if(files != None):
		pheno['files']=files
	pheno['meta_data']=meta_data

	return Phenopacket(**pheno)



def convert_time(time_datetime):
#    dt = datetime.strptime(episode_data['timestamp'], '%Y-%m-%d %H:%M:%S.%f')
	return datetime.strftime(time_datetime, '%Y-%m-%dT%H:%M:%S.%fZ')

def convert_time_toseconds(time_string):
	dt=datetime.strptime(time_string,'%Y-%m-%dT%H:%M:%S.%fZ')
	return int(datetime.timestamp(dt))


def main():

	parser = argparse.ArgumentParser()
	parser.add_argument('--loglevel',help='the logging level:DEBUG,INFO,WARNING,ERROR or CRITICAL',default='WARNING')
	parser.add_argument('--dbname',help='name of the db with omop cdr data',default='synthea10')
	parser.add_argument('--username',help='username for omop cdr db',default='postgres')  
	parser.add_argument('--password',help='password for omop cdr db',default='lollipop')
	parser.add_argument('--host',help='db host',default='172.23.0.2')
	parser.add_argument('--port',help='db port',default=5432)
   # parser.add_argument('--check',action='store_true', help='check the missing leafs for leafs that should be there but are not')
	parser.add_argument('--patient_id',help='id of the patient',default=1)
	parser.add_argument('--myname',help='phenotype creator name',default='creator')
	parser.add_argument('--outputfile',help='phenopacket output json file',default='pheno.output')

#    parser.add_argument('--outputfilebasename',help='output file basename',default='output')
	args=parser.parse_args()


	loglevel=getattr(logging, args.loglevel.upper(),logging.WARNING)
	if not isinstance(loglevel, int):
		raise ValueError('Invalid log level: %s' % loglevel)
	logging.basicConfig(filename='./conversion.log',filemode='w',level=loglevel)

	patient_id=args.patient_id
	myname=args.myname
	outputfile=args.outputfile

	# Connect to db
	dbname=args.dbname
	username=args.username
	password=args.password
	host=args.host
	port=args.port
	print(f'Connecting to db={dbname} host={host} port={port} with username={username} password={password}')
	print(f'Chosen Patient={patient_id}')
	logging.info(f'Connecting to db={dbname} host={host} port={port} with username={username} password={password}')
	logging.info(f'Chosen Patient={patient_id}')

	conn = get_postgres_connection(dbname,username,password,host,port)


	# Get Individual records
	cur = conn.cursor()
	cur.execute(get_individual_query(patient_id))
	records = cur.fetchall()


	logging.debug(f'Individual:records for patient {patient_id} =\n {records}')


	mydict=parse_Individual(records)

	logging.debug(f'mydict {mydict}')

	vsdict={}
	if(mydict['vital_status']!=0):
		#Get Vital Status records
		cur = conn.cursor()
		cur.execute(get_vitalstatus_query(patient_id))
		records = cur.fetchall()

		logging.debug(f'Vital Status:records for patient {patient_id} =\n {records}')


		vsdict=parse_VitalStatus(records)

		logging.debug(f'Vital Status dict {vsdict}')

	#Arrange dict for Pheno

	individualdict=createDictIndividual(mydict,vsdict)

	logging.debug(f'Whole Individual dict {individualdict}')


	#Get Measurments records
	cur = conn.cursor()
	cur.execute(get_measurement_query(patient_id))
	records = cur.fetchall()

	logging.debug(f'Measurement:records for patient {patient_id} =\n {records}')

	md=parse_Measurements(records)

	logging.debug(f'Measurement dictsList: {md}')

	mlist=createListDictMeasurements(md)

	logging.debug(f'mlist: {mlist}')

	#Create Pheno

	individualpheno=createPhenoIndividual(individualdict)

	measurementpheno=createPhenoMeasurement(mlist)

	meta_data=createMetadata(myname)

	pheno=createPheno(myid=patient_id,subject=individualpheno,measurements=measurementpheno,meta_data=meta_data)

	#Write Pheno to outputfile
	ijson= MessageToJson(pheno)
	logging.debug(f'pheno json: {ijson}')
	with open(outputfile,'w') as of:
		of.write(ijson)
	print(f'phenopacket written in file {outputfile}')
	logging.info(f'phenopacket written in file {outputfile}')

if __name__ == '__main__':
	main()
