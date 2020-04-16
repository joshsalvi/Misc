
# coding: utf-8

# # Extract text of Admit Notes
# #### Start with patient MRN and the Admission Date.
# #### This requires 4 sequential sql pulls:
# 1) Find the PatientEncounterID for the hospital admission.
# 2) Find the set of NoteIDs associated with the ZPatientID (PatientEncounterID missing from Notes.
# 3) Find set of NoteIDs associated with ZPatientID AND date of admission.
# 4) Find the text associated with the NoteIDs.
# After get notes, isolate Psychiatric Admission Note and Admission H&P
#

# In[1]:


import base64
import os
import sqlalchemy
import getpass
import sys
import datetime
import pandas as pd
import re
import numpy as np
from scipy import stats

userid = getpass.getuser()
print(userid)
pswd = getpass.getpass('Provide password to connect:')

connection_str ="mssql+pymssql://PARTNERS\\" + str(userid) + ":" + pswd + "@PHSEDW.partners.org"
engine = sqlalchemy.create_engine(connection_str)
#does not establish DBAPI connection UNTIL engine.connect called


# Input MRN, Admit Date

# In[4]:


MRNFile=pd.read_csv('EPIC_GetNotesEDW.csv')
MRNFile.head()

MRN = MRNFile['MRN']
#MRN.to_csv('MRNString.csv', index=False)
#MRNAll = open('MRNString.csv', 'r')
#mrnNote = MRNAll.read()
#mrnNote = mrnNote.splitlines()

#AdmitDate = MRNFile['DateAdmit']
AdmitDate = '2018-08-24'
#AdmitDate.to_csv('DateAdmit.csv', index=False)
#admitDateNote = open('DateAdmit.csv', 'r')
#dateNote = admitDateNote.read()
#dateNote = dateNote.splitlines()
#Don't need any of above

#only do few patients at time
lengthMRN = len(MRN)
lengthAdmitDate = len(AdmitDate)

#test first patient of 1,277 patients
#MRNCurrList = [number for number in range(0,1)]
#dateCurrList = [number for number in range(0,1)]
##After ran script and obtained all notes, 1,255 left (excluded few who were past admits)
MRNCurrList = [number for number in range(308,309)]
#dateCurrList = [number for number in range(0,lengthAdmitDate)]

print(MRN[MRNCurrList])
#print(AdmitDate[dateCurrList])
print(AdmitDate)


# ## Since can't split for loop across multiple cells, will define steps:
# ### 1) Find ZPatient ID (MCL tables do not always have PatientEncounterID for the hospital admission, which would be ideal)
# ### 2) Two inner joins: use ZPatient ID to get Note ID's associated with patient (Note Table). Link NoteIDs to Date of Service (NoteText Table)
# ### 3) Each patient has many notes. Select only McLean (112) and date of admission (AdmitDate). Get rid of duplicate notes (with same NoteID). Each patient will still have multiple notes
# ### 4) Get notes corresponding to Note ID's. Search for CEC Admission Note

# In[5]:


#1) Find ZPatientID and PatientEncounterID

#set up the query;
for whichPt in MRNCurrList:
    sql_string = """
    SELECT
      ptenchosp.PatientEncounterID,
      id.PatientIdentityID,
      id.PatientID,
      ptenchosp.ExpectedAdmissionDTS
    FROM
      Epic.Encounter.PatientEncounterHospital_McLean AS ptenchosp
    INNER JOIN
      Epic.Patient.Identity_McLean AS id
    ON
      id.PatientID = ptenchosp.PatientID
    WHERE
      id.IdentityTypeID = 112
      AND id.PatientIdentityID = '{}'
      AND CAST(ptenchosp.HospitalAdmitDTS AS date) = '{}'
    ORDER BY
      ptenchosp.HospitalAdmitDTS
    """

    #run the query, inserting the parameters into the query
    with engine.connect() as cxn:
        currMRN = MRN[whichPt]
        currAdmitDate = AdmitDate
        #currAdmitDate = AdmitDate[whichPt]

        print(currMRN)
        print(currAdmitDate)
        PtEncounterID = pd.io.sql.read_sql(sql_string.format(currMRN, currAdmitDate), cxn)
        #print(PtEncounterID)

    #display a warning if there were multiple admissions; try taking this out
    if len(PtEncounterID) > 1:
        warn_string = 'Warning: More than one admission for {} on {}. Using most recent admission on that date.'
        #print(warn_string.format(MRN, AdmitDate))

    ZPatientID = PtEncounterID.iloc[0]['PatientID']
    print(ZPatientID)

    #pick out the PatientEncounterID
    PtEncounterID = PtEncounterID.iloc[0]['PatientEncounterID'] #Use index 0 for first admission; -1 for last admission
    PtEncounterID = int(PtEncounterID)
    #print(PtEncounterID)

    #2. Two inner joins: use ZPatient ID to get Note ID's associated with patient (Note Table).
    #Link NoteIDs to Date of Service (NoteText Table)
    #set up the query
    sql_string2 = """
    SELECT
        notes.NoteID,
        id.PatientID,
        id.PatientIdentityID,
        id.IdentityTypeID,
        notetext.ContactDTS
    FROM
        Epic.Clinical.Note_McLean AS notes
        INNER JOIN Epic.Patient.Identity_McLean AS id ON id.PatientID = notes.PatientLinkID
        INNER JOIN Epic.Clinical.NoteText_McLean AS notetext ON notes.NoteID = notetext.NoteID
    WHERE
        notes.PatientLinkID = '{}'
    ORDER BY
        notes.NoteID
    """

    #print(sql_string2)


    #run the query, inserting the parameters into the query
    with engine.connect() as cxn:
        NoteID = pd.io.sql.read_sql(sql_string2.format(ZPatientID), cxn)

    #found there were many duplicate NoteID's for some patients
    #3. Convert to dataframe.
    #Next use dataframe Note ID to select McLean notes (112) and date of admission (AdmitDate)
    #Get rid of duplicates (keep first)
    NoteIDFrame = pd.DataFrame(NoteID)
    #get rid of non-McLean notes first
    NoteIDFrame = NoteIDFrame.where(NoteIDFrame['IdentityTypeID'] == 112.0)
    NoteIDFrame['ContactDTS'] = pd.to_datetime(NoteIDFrame['ContactDTS'], unit='s')
    NoteIDFrame = NoteIDFrame.where(NoteIDFrame['ContactDTS'] == AdmitDate)
    NoteIDFrame = NoteIDFrame.dropna()
    NoteIDFrame = NoteIDFrame.drop_duplicates(subset='NoteID', keep='first')
    #sort by Note ID
    NoteIDFrame = NoteIDFrame.sort_values(by='NoteID')
    #renumber indices, drop=True gets rid of old indices
    NoteIDFrame = NoteIDFrame.reset_index(drop=True)
    #print(NoteIDFrame)

    #get list of note ID for patient
    listNoteID = list(NoteIDFrame['NoteID'])
    #determine number of notes for patient that occurred on Day of Admit
    numberNotes = len(listNoteID)
    #print(listNoteID)

    #4) Get Notes corresponding to Note ID's
    #set up the query
    sql_string = """
    SELECT
      NoteTXT
    FROM
      Epic.Clinical.NoteText_McLean
    WHERE
      NoteID = '{}'
    ORDER BY
      LineNBR
    """
    #print(sql_string)

    #run the query, inserting the parameters into the query
    #filename MRN_NoteID
    #search each note for Medical Admission Note and Psychiatric Admission Note in first line
    noteCounter = 0
    for patientList in listNoteID:
        if noteCounter < 6:
            with engine.connect() as cxn:
                NoteText = pd.io.sql.read_sql(sql_string.format(patientList), cxn)
                fulltext = NoteText.NoteTXT.str.cat()
                filename = [str(MRN[whichPt]) +'_' + str(patientList) +'.txt']
                filename = "".join(filename)
                #print(filename)
                f = open(filename, 'w')
                f.write(fulltext)
                f.close()
                f = open(filename, 'r')
                CECnote = f.readline()
                psychNote = re.findall('McLean Clinical Evaluation Center  Psychiatric Admission Note', CECnote)
                medNote = re.findall('CEC Medical Admission Note', CECnote)
                if len(psychNote) > 0:
                    noteCounter = noteCounter + 1
                    psychFileName = ['PsychAdmit_' + str(MRN[whichPt]) + '.txt']
                    psychFileName = "".join(psychFileName)
                    print(psychFileName)
                    os.rename(filename, psychFileName)
                    #f = open(psychFileName, 'w')
                    #f.write(fulltext)
                    #f.close()
                if len(medNote) > 0:
                    noteCounter = noteCounter + 1
                    medFileName = ['MedAdmit_' +str(MRN[whichPt]) +'.txt']
                    medFileName = "".join(medFileName)
                    print(medFileName)
                    os.rename(filename, medFileName)
                    #f = open(medFileName, 'w')
                    #f.write(fulltext)
                    #f.close()
