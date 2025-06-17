"""
Predict with the HiDRA model

Requirement:
    model.hdf: Pre-trained model file
    predict.csv: The list of cancer-drug pairs that will be predicted. Consists of idx (int), Drug name, Cnacer cell line name.
    input_dir: The directory that includes input files
"""
#Import basic packages
import numpy as np
import pandas as pd
import csv

import os
import argparse

#Import keras modules
import tensorflow as tf
import keras.backend as K
import keras.backend.tensorflow_backend as KTF
import keras
import keras.layers
from keras.layers import Layer 
import keras.initializers
from keras.models import Model, Sequential,load_model
from keras.layers import Input, Dense, Dropout, BatchNormalization, Activation, Multiply, multiply,dot
from keras.layers import Concatenate,concatenate
from keras.optimizers import Adam
from keras.utils import plot_model

import re
from utils_def import *
from utils_hidra import *

# #Fix the random seed
# np.random.seed(5)

import pickle
with open("_data/geneset.pickle", 'rb') as f :
    GeneSet_Dic_withoutNA = pickle.load(f)


#Using 20% of GPU memory only
def get_session(gpu_fraction=0.1):

    num_threads = os.environ.get('OMP_NUM_THREADS')
    gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=gpu_fraction)

    if num_threads:
        return tf.Session(config=tf.ConfigProto(gpu_options=gpu_options, intra_op_parallelism_threads=num_threads))
    else:
        return tf.Session(config=tf.ConfigProto(gpu_options=gpu_options))

# def read_files(Predict_list,Input_directory):
#     
#     #Read predict list
#     GDSC_without_duplicated=pd.read_csv(Predict_list,index_col=0)
#     GDSC_without_duplicated.columns=['Drug name','Cell line name']
#     
#     #Read input files for HiDRA model
#     X_origin=[pd.read_csv(Input_directory+str(x)+'.csv',index_col=0) for x in range(186)]
#     X_origin=[df.loc[GDSC_without_duplicated['Cell line name']] for df in X_origin]
#     Drug=pd.read_csv(Input_directory+'drug.csv',index_col=0)
#     X_origin.append(Drug.loc[GDSC_without_duplicated['Drug name']])
#     
#     return X_origin


def read_files(Input_directory, file_drug="drug.csv"):

    #Read input files for HiDRA model
    X_origin = [pd.read_csv(Input_directory+"/"+str(x)+'.csv',index_col=0) for x in range(len(GeneSet_Dic_withoutNA))]
    Drug = pd.read_csv(Input_directory+"/"+file_drug,index_col=0)
    X_origin.append(Drug) 

    return X_origin


#Make new HiDRA model
def Making_Model():
    # #Read Gene expression and Gene set for model making
    # #They are same with the source codes in the function 'Read_files'
    # #Read Gene expression file
    # GeneExpression_with_Symbol=pd.read_csv("expression.csv")
    # GeneExpression_with_Symbol.index=GeneExpression_with_Symbol.Gene_Symbol
    # 
    # 
    # #Read Gene set file
    # GeneSet_List=[]
    # with open("geneset.gmt") as f:
    #     reader = csv.reader(f)
    #     data = list(list(rec) for rec in csv.reader(f, delimiter='\t')) #reads csv into a list of lists
    #     for row in data:
    #         GeneSet_List.append(row)
    # 
    # GeneSet_Dic={}
    # for GeneSet in GeneSet_List:
    #     GeneSet_Dic[GeneSet[0]]=GeneSet[2:]
    # 
    # GeneSet_Dic_withoutNA={}
    # for GeneSet in GeneSet_Dic:
    #     GeneSet_Dic_withoutNA[GeneSet]=GeneExpression_with_Symbol['ES3'][GeneSet_Dic[GeneSet]].dropna().index.values

    #HiDRA model with keras
    #Drug-level network
    Drug_feature_length=512
    Drug_Input=Input((Drug_feature_length,), dtype='float32', name='Drug_Input')
    
    Drug_Dense1=Dense(256, name='Drug_Dense_1')(Drug_Input)
    Drug_Dense1=BatchNormalization(name='Drug_Batch_1')(Drug_Dense1)
    Drug_Dense1=Activation('relu', name='Drug_RELU_1')(Drug_Dense1)

    Drug_Dense2=Dense(128, name='Drug_Dense_2')(Drug_Dense1)
    Drug_Dense2=BatchNormalization(name='Drug_Batch_2')(Drug_Dense2)
    Drug_Dense2=Activation('relu', name='Drug_RELU_2')(Drug_Dense2)
    
    #Drug network that will be used to attention network in the Gene-level network and Pathway-level network
    Drug_Dense_New1=Dense(128, name='Drug_Dense_New1')(Drug_Input)
    Drug_Dense_New1=BatchNormalization(name='Drug_Batch_New1')(Drug_Dense_New1)
    Drug_Dense_New1=Activation('relu', name='Drug_RELU_New1')(Drug_Dense_New1)

    Drug_Dense_New2=Dense(32, name='Drug_Dense_New2')(Drug_Dense_New1)
    Drug_Dense_New2=BatchNormalization(name='Drug_Batch_New2')(Drug_Dense_New2)
    Drug_Dense_New2=Activation('relu', name='Drug_RELU_New2')(Drug_Dense_New2)

    #Gene-level network
    GeneSet_Model=[]
    GeneSet_Input=[]
    
    #Making networks whose number of node is same with the number of member gene in each pathway    
    for GeneSet in GeneSet_Dic_withoutNA.keys():
        Gene_Input=Input(shape=(len(GeneSet_Dic_withoutNA[GeneSet]),),dtype='float32', name=GeneSet+'_Input')
        Drug_effected_Model_for_Attention=[Gene_Input]
        #Drug also affects to the Gene-level network attention mechanism
        Drug_Dense_Geneset=Dense(int(len(GeneSet_Dic_withoutNA[GeneSet])/4)+1,dtype='float32',name=GeneSet+'_Drug')(Drug_Dense_New2)
        Drug_Dense_Geneset=BatchNormalization(name=GeneSet+'_Drug_Batch')(Drug_Dense_Geneset)
        Drug_Dense_Geneset=Activation('relu', name=GeneSet+'Drug_RELU')(Drug_Dense_Geneset)
        Drug_effected_Model_for_Attention.append(Drug_Dense_Geneset) #Drug feature to attention layer
 
        Gene_Concat=concatenate(Drug_effected_Model_for_Attention,axis=1,name=GeneSet+'_Concat')
        #Gene-level attention network
        Gene_Attention = Dense(len(GeneSet_Dic_withoutNA[GeneSet]), activation='tanh', name=GeneSet+'_Attention_Dense')(Gene_Concat)
        Gene_Attention=Activation(activation='softmax', name=GeneSet+'_Attention_Softmax')(Gene_Attention)
        Attention_Dot=dot([Gene_Input,Gene_Attention],axes=1,name=GeneSet+'_Dot')
        Attention_Dot=BatchNormalization(name=GeneSet+'_BatchNormalized')(Attention_Dot)
        Attention_Dot=Activation('relu',name=GeneSet+'_RELU')(Attention_Dot)
        
	      #Append the list of Gene-level network (attach new pathway)
        GeneSet_Model.append(Attention_Dot)
        GeneSet_Input.append(Gene_Input)

    Drug_effected_Model_for_Attention=GeneSet_Model.copy()
    
    #Pathway-level network
    Drug_Dense_Sample=Dense(int(len(GeneSet_Dic_withoutNA)/16)+1,dtype='float32',name='Sample_Drug_Dense')(Drug_Dense_New2)
    Drug_Dense_Sample=BatchNormalization(name=GeneSet+'Sample_Drug_Batch')(Drug_Dense_Sample)
    Drug_Dense_Sample=Activation('relu', name='Sample_Drug_ReLU')(Drug_Dense_Sample)    #Drug feature to attention layer
    Drug_effected_Model_for_Attention.append(Drug_Dense_Sample)
    GeneSet_Concat=concatenate(GeneSet_Model,axis=1, name='GeneSet_Concatenate')
    Drug_effected_Concat=concatenate(Drug_effected_Model_for_Attention,axis=1, name='Drug_effected_Concatenate')
    #Pathway-level attention
    Sample_Attention=Dense(len(GeneSet_Dic_withoutNA.keys()),activation='tanh', name='Sample_Attention_Dense')(Drug_effected_Concat)
    Sample_Attention=Activation(activation='softmax', name='Sample_Attention_Softmax')(Sample_Attention)
    Sample_Multiplied=multiply([GeneSet_Concat,Sample_Attention], name='Sample_Attention_Multiplied')
    Sample_Multiplied=BatchNormalization(name='Sample_Attention_BatchNormalized')(Sample_Multiplied)
    Sample_Multiplied=Activation('relu',name='Sample_Attention_Relu')(Sample_Multiplied)
    
    #Making input list
    Input_for_model=[]
    for GeneSet_f in GeneSet_Input:
        Input_for_model.append(GeneSet_f)
    Input_for_model.append(Drug_Input)
    
    #Concatenate two networks: Pathway-level network, Drug-level network 
    Total_model=[Sample_Multiplied,Drug_Dense2]
    Model_Concat=concatenate(Total_model,axis=1, name='Total_Concatenate')

    #Response prediction network
    Concated=Dense(128, name='Total_Dense')(Model_Concat)
    Concated=BatchNormalization(name='Total_BatchNormalized')(Concated)
    Concated=Activation(activation='relu', name='Total_RELU')(Concated)

    Final=Dense(1, name='Output')(Concated)
    model=Model(inputs=Input_for_model,outputs=Final)
    
    return model



def main():
    KTF.set_session(get_session())

    #Reading argument 
    parser=argparse.ArgumentParser(description='HiDRA:Hierarchical Network for Drug Response Prediction with Attention-Predict')
    
    # #Options
    # parser.add_argument('-m',type=str,help='The model file')
    # parser.add_argument('-p',type=str,help='The prediction list')
    # parser.add_argument('-i',type=str,help='The input directory')
    # parser.add_argument('-o',type=str,help='The output file path that prediction result be stored')
    
    parser.add_argument('-ic50', type=str, default="IC50_GDSC2.txt", help='IC50 Files')
    parser.add_argument('-dir_in', type=str, default="_data/ProcessedFile", help='Input directory')
    parser.add_argument('-drug', type=str, default="drug_chembl.csv", help='Input directory')
    
    parser.add_argument("-cpu", type=int, default=4)
    # parser.add_argument('-gpu', type=int, default=0, help='Use GPU/CPU')

    parser.add_argument("-col_cell", type=str, default="COSMIC_ID")
    parser.add_argument("-col_drug", type=str, default="Drug")
    parser.add_argument("-col_ic50", type=str, default="LN_IC50")
    
    parser.add_argument("-dir_test", type=str, default=None)
    parser.add_argument("-dir_param", type=str, default=None)
    
    args=parser.parse_args()
    
    # gpu = args.gpu
    ic50 = args.ic50

    dir_ic50 = "./_data"
    sep = "\t" if ".txt" in args.ic50 else ","
    input_ic50 = os.path.join(dir_ic50, args.ic50)
    ic50_data = pd.read_csv(input_ic50, header=0, sep=sep)
    
    print("\n### Data Processing")
    data = read_files(args.dir_in, args.drug)
    data = index_to_str(data)
    ic50_data = filt_ic50_(ic50_data, data, args)
    
    batch_size = 1024
    test_loader = DataLoader(ic50_data, data, args, batch_size=batch_size, shuffle=False)
    
    model = Making_Model()
    model.load_weights(args.dir_param)
    
    
    # Test Performance
    print("\n### Test Performance...")
    import time
    start_time = time.perf_counter()
    pred_test = model.predict_generator(test_loader, workers=args.cpu)
    end_time = time.perf_counter()
    
    test_time = 1000*(end_time - start_time)
    # pred_to_csv(pred_test, ic50_data, args.dir_test) 
    
    
    # Inference Time
    if "ChEMBL" in args.ic50 and args.col_cell=="COSMIC_ID":
        args.seed = args.dir_test.split("seed")[-1].split(".")[0]
        args.dir_out = "/".join(args.dir_test.split("/")[:-1])
        args.dir_time = "{}/log_time_chembl_seed{}.csv".format(args.dir_out, args.seed)
        time_to_csv(test_time=test_time, dir_time=args.dir_time)
    
    # #Read input files from predict list 
    # model_input=read_files(args.p,args.i)
    # #Load model in hdf5 file format
    # model=load_model(args.m,compile=False)
    # #Predict
    # result=model.predict(model_input)
    # result=[y[0] for y in result]
    # predict_list=pd.read_csv(args.p)
    # predict_list['result']=result
    # #Save the predict results to the output directory
    # predict_list.to_csv(args.o)
    
if __name__=="__main__":
    main()
