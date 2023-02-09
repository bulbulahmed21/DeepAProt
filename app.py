##################################
############### Requirement----- tensroflow version 2.3.0...--------##########################
###################################

from _datetime import date, datetime
import json
import xlsxwriter
import os
import numpy as np
import pandas as pd
import tensorflow as tf
from keras.models import load_model, Sequential
from sklearn.preprocessing import StandardScaler
from flask import Flask, render_template, request, send_file
from werkzeug.utils import secure_filename
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv
from io import StringIO


# Paths to folders
upload_filepath_csv= os.path.join('custom_directory','uploaded_files','csv_files')
upload_filepath_fasta= os.path.join('custom_directory','uploaded_files','fasta_files')
model_filepath = os.path.join('custom_directory','saved_models')
model_results = os.path.join('custom_directory','results')

# extension allowed
ALLOWED_EXTENSIONS = {'fasta'}

# Flask App
app = Flask(__name__)

app.config['upload_filepath_csv'] = upload_filepath_csv
app.config['upload_filepath_fasta'] = upload_filepath_fasta
app.config['model_filepath'] = model_filepath
app.config['model_results'] = model_results

# functions
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# activation functions
def sielu(x):
    return 0.5 * x * (2 * tf.sigmoid(2 * tf.sqrt(2 / np.pi) * (x + 0.044715 * tf.pow(x, 3))))
    #return 0.5 * x * (1 + tf.sigmoid(tf.sqrt(2 / np.pi) * (x + 0.044715 * tf.pow(x, 3))))

def gelu(x):
    return 0.5 * x * (1 + tf.tanh(tf.sqrt(2 / np.pi) * (x + 0.044715 * tf.pow(x, 3))))

# Dataset Prepocessing
def data_preprocess(f):
    filename = f
    data = pd.read_csv(os.path.join(app.config['upload_filepath_csv'], filename), header=None)
    # print(data)
    # corr = data.corr()
    # columns = np.full((corr.shape[0],), True, dtype=bool)
    # for i in range(corr.shape[0]):
    #     for j in range(i + 1, corr.shape[0]):
    #         if corr.iloc[i, j] >= 0.8:
    #             if columns[j]:
    #                 columns[j] = False
    # selected_columns = data.columns[columns]
    # data_predict = data[selected_columns]
    seq_id = data.iloc[:, 0].values
    # print(seq_id)
    x = data.iloc[:, 1:].values  # all rows and 2nd to all remaining coloumn

    preprocess = StandardScaler()
    x = preprocess.fit_transform(x)
    x = np.reshape(x,(x.shape[0], x.shape[1], 1))
    return x, seq_id

#### function that extract features from fasta file to a csv file
def fastaToCSV(f):
    fasta_file=f
    f, e = os.path.splitext(fasta_file)
    fasta_file_path=os.path.join(app.config['upload_filepath_fasta'], fasta_file)
    # print(fasta_file)
    complete_feature_List = []
    # print(SeqIO.parse(fasta_file_path, "fasta"))
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        ### print(record)
        feature_row_list = []
        ##### Amino Acids
        seqid = str(record.id)
        # print(seqid)
        feature_row_list.append(seqid)
        seq = str(record.seq)
        X = ProteinAnalysis(seq)
        aaDict = X.count_amino_acids()
        for aa in aaDict:
            feature_row_list.append(aaDict.get(aa, 0))

        # percentage
        PDict = X.get_amino_acids_percent()
        for aa in PDict:
            feature_row_list.append(PDict.get(aa, 0))

        # Molecular Weight
        W = X.molecular_weight()
        feature_row_list.append(W)

        # Aromaticity
        A = X.aromaticity()
        feature_row_list.append(A)

        # Instability
        I = X.instability_index()
        feature_row_list.append(I)

        # GRAVY
        G = X.gravy()
        feature_row_list.append(G)

        # ISO-Electric point
        PI = X.isoelectric_point()
        feature_row_list.append(PI)

        # secondary
        SS = X.secondary_structure_fraction()
        feature_row_list.append(SS[0])
        feature_row_list.append(SS[1])
        feature_row_list.append(SS[2])

        complete_feature_List.append(feature_row_list)
    csv_file = f+'.csv'
    with open(os.path.join(app.config['upload_filepath_csv'], csv_file), 'w') as f11:
        writer = csv.writer(f11)
        writer.writerows(complete_feature_List)
    return csv_file


# Home page
@app.route('/',methods = ['POST', 'GET'])
def home_screen():
    return render_template('test.html')

# Dataset upload
@app.route('/upload',methods = ['POST', 'GET'])
def upload_file():
    if request.method == 'POST':
        # check if the post request has the file part
        #print(request.files)
        #print(request.form)
        if  request.form and request.files is None:
            #print('1')
            # flash('No file part')
            return json.dumps({
                'status':'NOT OK'
                               })
        elif 'file' in request.files:
            # print('1')
            datafile = request.files['file']
            # print(datafile)
            if (allowed_file(datafile.filename)):
                # print('2')
                if datafile:
                    # print('3')
                    filename = secure_filename(datafile.filename)
                    f, e = os.path.splitext(filename)
                    today = datetime.now()
                    today = today.strftime("%d_%m_%y_%H_%M_%S")
                    newFileNme = f + '_' + today + e
                    datafile.save(os.path.join(app.config['upload_filepath_fasta'], newFileNme))
                    new_csv_file = fastaToCSV(newFileNme)
                    return json.dumps({
                        'status': 'OK',
                        'file_name': new_csv_file
                    })
                else:
                    return json.dumps({
                        'status': 'NOT OK'
                    })
            else:
                return json.dumps({
                    'status': 'NOT CSV'
                })
        elif 'sequence' in request.form:
            #print('in sequence')
            seq = request.form['sequence']
            seq_str = StringIO(seq)
            #print(type(seq))
            #seq_df = pd.DataFrame(seq)
            today = datetime.now()
            today = today.strftime("%d_%m_%y_%H_%M_%S")
            seq_filename = 'pasted_seq_' + today + '.fasta'

            file = open(os.path.join(app.config['upload_filepath_fasta'], seq_filename), "w+")
            # out = '\n'.join(['>sp' + str(i + 1) + "\n" + j for i, j in enumerate(seq)])
            file.write(seq)
            file.close
            # seq.save(os.path.join(app.config['upload_filepath_fasta'], seq_filename))
            new_csv_file = fastaToCSV(seq_filename)
            return json.dumps({
                'status': 'OK',
                'file_name': new_csv_file
            })
        else:

            return json.dumps({
                'status': 'NOT OK'
            })
    else:
        return json.dumps({
            'status': 'NOT OK'
        })

# Prediction with the saved models
@app.route('/model',methods = ['POST', 'GET'])
def loading_model():
    filename = None
    # modeltype= None
    # activationtype = None
    # act=None
    # mod=None
    if request.method == 'POST':
        # print("okk")
        if  request.form is None:
            # print('in_model_none')
            # flash('No file part')
            return json.dumps({
                'status':'NOT OK'
                               })
        else:
            # print('in_model_file')
            filename = request.form['file_name']
            # print(filename)
            f, e = os.path.splitext(filename)
            x, seq_id = data_preprocess(filename)

            final_result = pd.DataFrame(seq_id, columns=['Sequence_ID'])
            # print(final_result.head())

            my_model_1 = load_model(os.path.join(app.config['model_filepath'], 'model_cold_sielu.h5'),
                                    custom_objects={'sielu': sielu})
            my_model_2 = load_model(os.path.join(app.config['model_filepath'], 'model_drought_sielu.h5'),
                                                             custom_objects={'sielu': sielu})
            my_model_3 = load_model(os.path.join(app.config['model_filepath'], 'model_heat_sielu.h5'),
                                                             custom_objects={'sielu': sielu})
            my_model_4 = load_model(os.path.join(app.config['model_filepath'], 'model_salt_sielu.h5'),
                                                             custom_objects={'sielu': sielu})
            models =(my_model_1,my_model_2, my_model_3, my_model_4)
            col_list = ['Cold Stress', 'Drought Stress', 'Heat Stress', 'Salt Stress']
            for m in models:
                y_pred = m.predict_classes(x)
                # print(type(y_pred))
                new_y = []
                for i in y_pred:
                    # print(y_pred[i])
                    if (y_pred[i] == 1):
                        new_y.append('Present')
                    else:
                        new_y.append('Absent')
                # print(new_y)
                if(m==my_model_1):
                    final_result['Cold_Stress'] = new_y
                    # print(final_result.head())
                elif(m==my_model_2):
                    final_result['Drought_Stress'] = new_y
                    # print(final_result.head())
                elif(m==my_model_3):
                    final_result['Heat_Stress'] = new_y
                    # print(final_result.head())
                elif(m==my_model_4):
                    final_result['Salt_Stress'] = new_y
                    # print(final_result.head())
                else:
                    print()
            print(final_result)
            data = final_result.to_json(orient='records')
            # print(data)
            return data
    else:
        return json.dumps({
            'status': 'NOT OK'
        })

# File download
@app.route('/download/<file>',methods = ['GET'])
def download(file):
    path = os.getcwd()
    path_rs=path+"/"+app.config['model_results']+ "/"+file
    return send_file(path_rs,
                    mimetype='text/csv',
                    attachment_filename=file,
                    as_attachment=True)

# run the api
if __name__ == '__main__':
   app.run(host= "172.16.10.3", port=5000, threaded= False)
