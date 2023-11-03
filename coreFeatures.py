########################################################################################
import os
from os import path as p
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
import importlib
import subprocess
import argpass
from custom_utils import pdb2df, area2df
########################################################################################
########################################################################################
# get inputs
def read_inputs():
    parser = argpass.ArgumentParser()
    parser.add_argument("--config")
    args = parser.parse_args()
    configName=args.config
    if  args.config == None:
        print('No config file name provided.')
    configName = p.splitext(configName)[0]

    # add config to PYTHONPATH
    cwd = os.getcwd()
    configPath = p.join(cwd,configName)
    sys.path.append(configPath)
    # import config file and run input function to return variables
    try:
        config_module = __import__(configName)
        inputDir, outputDir, msmsDir, msmsExe, pdb2xyzrExe, aminoAcidTable = config_module.inputs()
        return inputDir, outputDir, msmsDir, msmsExe, pdb2xyzrExe, aminoAcidTable
    except ImportError:
        print(f"Error: Can't to import module '{configName}'. Make sure the input exists!")
        print("HOPE IS THE FIRST STEP ON THE ROAD TO DISAPPOINTMENT")
        exit()

########################################################################################

def initialiseAminoAcidInformation(aminoAcidTable):
    AminoAcidNames = ["ALA","ARG","ASN","ASP","CYS",
                      "GLN","GLU","GLY","HIS","ILE",
                      "LEU","LYS","MET","PHE","PRO",
                      "SER","THR","TRP","TYR","VAL"]  

    # read file with amino acids features
    aminoAcidProperties = pd.read_csv(
        aminoAcidTable, sep="\t", index_col=1
    )
    aminoAcidProperties.index = [el.upper() for el in aminoAcidProperties.index]
    aminoAcidProperties = aminoAcidProperties.iloc[:, 1:]

    return AminoAcidNames, aminoAcidProperties

########################################################################################
def getPdbList(dir):
    pdbList=[]
    idList=[]
    for file in os.listdir(dir):
        fileData = p.splitext(file)
        if fileData[1] == '.pdb':
            idList.append(fileData[0])
            pdbList.append(p.join(dir,file))
    return idList, pdbList

########################################################################################
def findCoreExterior(pdbFile,msmsDir,pdbDf,proteinName,outDir):
    # change working directory so MSMS can find all the files it needs
    os.chdir(msmsDir)
    # find executables
    pdb2xyzrExe = "./pdb_to_xyzr"
    msmsExe = "./msms.x86_64Linux2.2.6.1"
    # convert pdb file to MSMS xyzr file 
    xyzrFile = p.join(outDir, f'{proteinName}.xyzr')
    command = f"{pdb2xyzrExe} {pdbFile} > {xyzrFile}"
    subprocess.run(command, shell=True)
    # use MSMS to create an area file
    areaOut = p.join(outDir,proteinName)
    command = f"{msmsExe} -if {xyzrFile} -af {areaOut}"
    subprocess.run(command, shell=True)
    areaFile=p.join(outDir,f"{proteinName}.area")
    # convert area file to dataframe, merge with main pdb dataframe
    areaDf = area2df(areaFile=areaFile)
    pdbDf = pd.concat([pdbDf,areaDf],axis=1)

    # Group by residue and calculate the average SES score
    meanSesPerResidue = pdbDf.groupby('RES_SEQ')['SES'].mean()

    # Get residue sequences with average SES > 1
    exteriorResiduesIndex = meanSesPerResidue[meanSesPerResidue > 1].index

    # Split the DataFrame based on average SES > 1
    exteriorDf = pdbDf[pdbDf['RES_SEQ'].isin(exteriorResiduesIndex)]
    coreDf = pdbDf[~pdbDf['RES_SEQ'].isin(exteriorResiduesIndex)]

    # clean up
    os.remove(xyzrFile)
    os.remove(areaFile)

    return exteriorDf, coreDf
########################################################################################
def get_counts_in_region(coreDf,extDf,pdbDf,proteinName,aminoAcidNames):
    featuresDict={}
    for region, df in zip(["core","ext","protein"],[coreDf,extDf,pdbDf]):
        # initialise a features dataframe with column names
        columnNames=[]
        columnNames.append(f'{region}.total')
        for aminoAcid in aminoAcidNames:
            columnNames.append(f'{region}.{aminoAcid}')
        for element in ["C","N","O"]:
            columnNames.append(f'{region}.{element}')
        columnNames=columnNames.sort()
        countsDf = pd.DataFrame(columns=columnNames,index=[proteinName])

        # count elements [Carbon, Nitrogen, Oxygen]
        for element in ["C","N","O"]:
            try:
                countsDf.loc[:,f'{region}.{element}'] = df["ELEMENT"].value_counts()[element]
            except:
                countsDf.loc[:,f'{region}.{element}'] = 0
        # get unique residues only
        uniqueResDf= df.drop_duplicates(subset=["RES_SEQ"])
        totalRes=0
        # add get residue counts
        for aminoAcid in aminoAcidNames:
            if aminoAcid in df["RES_NAME"].unique():
                countsDf.loc[:,f'{region}.{aminoAcid}'] = uniqueResDf["RES_NAME"].value_counts()[aminoAcid]
            else:
                countsDf.loc[:,f'{region}.{aminoAcid}'] = 0
            totalRes+=countsDf[f'{region}.{aminoAcid}']
        countsDf[f"{region}.total"] = totalRes
        featuresDict.update({f'{region}.counts':countsDf})
    return featuresDict
########################################################################################
def get_properties_in_region(featuresDict, proteinName, aminoAcidNames, aminoAcidProperties):
    # loop through three regions
    for region in ["core","ext","protein"]:
        # initialise dataframe with columnNames
        columnNames=[]
        for property in aminoAcidProperties.columns:
            columnNames.append(f'{region}.{property}')
        propertiesDf=pd.DataFrame(columns=columnNames,index=[proteinName])
        # find count dataframe for corresponding region
        countDf = featuresDict[f'{region}.counts']
        for property in aminoAcidProperties:
            propertyValue=0
            for aminoAcid in aminoAcidNames:
                aaCount = countDf.at[proteinName,f"{region}.{aminoAcid}"]
                aaPropertyvalue = aminoAcidProperties.at[aminoAcid,property]
                value = aaCount * aaPropertyvalue
                propertyValue += value 

            totalAminoAcids=countDf.at[proteinName,f'{region}.total']
            if not totalAminoAcids == 0:
                propertyValue = propertyValue / totalAminoAcids
            propertiesDf[f'{region}.{property}'] = propertyValue
        featuresDict.update({f'{region}.properties':propertiesDf})
    return featuresDict
########################################################################################
def process_pdbs(pdbList, outDir, aminoAcidNames, aminoAcidProperties,msmsDir):
    allFeaturesList=[]
    for pdbFile in pdbList:
        proteinName=p.splitext(p.basename(pdbFile))[0]
        pdbDf = pdb2df(pdbFile=pdbFile)
        exteriorDf, coreDf = findCoreExterior(pdbFile=pdbFile, pdbDf=pdbDf,
                         proteinName=proteinName, msmsDir=msmsDir,
                         outDir=outDir)
        featuresDict = get_counts_in_region(coreDf=coreDf,extDf=exteriorDf,pdbDf=pdbDf,
                                                                    proteinName=proteinName,
                                                                    aminoAcidNames=aminoAcidNames)
        featuresDict = get_properties_in_region(featuresDict=featuresDict,
                                                proteinName=proteinName,
                                                aminoAcidNames=aminoAcidNames,
                                                aminoAcidProperties=aminoAcidProperties)
        featuresDf = pd.concat(featuresDict.values(),axis=1)
        allFeaturesList.append(featuresDf)
    allFeaturesDf=pd.concat(allFeaturesList,axis=0)

    allFeaturesDf.to_csv(p.join(outDir,"core_features.csv"),index=True)

########################################################################################
def main():
    # load user inputs
    read_inputs()
    outDir=os.makedirs(outputDir, exist_ok=True)
    # initialise amino acid data
    aminoAcidNames, aminoAcidProperties = initialiseAminoAcidInformation(aminoAcidTable)
    # get list of pdbFiles in pdbDir
    idList, pdbList = getPdbList(inputDir)

    # loop through permutations of orb and {region}.values
    process_pdbs(pdbList=pdbList,
                    outDir=outDir, 
                    aminoAcidNames=aminoAcidNames, aminoAcidProperties=aminoAcidProperties,
                    msmsDir=msmsDir)

    print("\nAll features have been generated and saved!")                
########################################################################################
main()