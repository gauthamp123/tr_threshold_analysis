import pandas as pd
import numpy as np
import math
from pathlib import Path

import pandas as pd
from pyfaidx import Fasta
import random

import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Conv1D, MaxPooling1D, GlobalMaxPooling1D, Dense, Concatenate, Dropout
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

tf.keras.mixed_precision.set_global_policy('mixed_float16')

def prepare_motif_encodings(sequence):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    encoded = np.zeros((len(sequence), 4))
    for i, base in enumerate(seq[:len(sequence)]):
        if base in mapping:
            encoded[i, mapping[base]] = 1
    return encoded



if __name__ == "__main__":
    df = pd.read_csv('data/expression_tr_annotations.txt', sep='\t')
    df.drop('Unnamed: 0', axis=1, inplace=True)


    columns = df.columns
    gene_exp_cols = [x for x in columns if 'GTEX' in x and ('Cortex' in x or 'cortex' in x)]
    numeric_features = [y for y in columns if 'Brain' in y and y not in gene_exp_cols and 'variant' not in y and ('Cortex' in y or 'cortex' in y)]
    max_len = (df['ann_end'] - df['ann_start']).max()
    
    fasta = Fasta('data/hg38.ml.fa')
    sequences = []
    disease = []
    max_len = (df['ann_end'] - df['ann_start']).max()
    sequence_encodings = []
    for idx, row in df.iterrows():
        start, end, chr = row['ann_start'], row['ann_end'], row['chr']
        seq = str(fasta[chr][start - 1 : end])
        sequences.append(seq)
        sequence_encodings.append(prepare_motif_encodings(seq))
        if isinstance(row['disease'], str):
            disease.append(1)
        else:
            disease.append(0)
        
    
    

    df['sequence'] = sequences
    df['disease_binary'] = disease

    breakpoint()
    scaler = StandardScaler()
    df[numeric_features] = scaler.fit_transform(df[numeric_features])

    X_seq = np.array(sequence_encodings)
    X_numeric = df[numeric_features].values
    y = np.array(df['disease_binary'])
    X_seq_train, X_seq_test, X_numeric_train, X_numeric_test, y_train, y_test = train_test_split(
        X_seq, X_numeric, y, test_size=0.2, random_state=42
    )

    seq_input = Input(shape=(None, 4), name="sequence_input")
    x = Conv1D(filters=32, kernel_size=10, activation="relu")(seq_input)
    x = GlobalMaxPooling1D()(x)  # Handles variable-length output


    # Numeric feature input
    numeric_input = Input(shape=(len(numeric_features),), name="numeric_input")
    y = Dense(64, activation="relu")(numeric_input)
    y = Dropout(0.3)(y)

    combined = Concatenate()([x, y])
    z = Dense(128, activation="relu")(combined)
    z = Dropout(0.4)(z)
    output = Dense(1, activation="sigmoid", name="output")  # Binary classification

    model = Model(inputs=[seq_input, numeric_input], outputs=output)
    model.compile(optimizer="adam", loss="binary_crossentropy", metrics=["accuracy"])
    model.summary()

    print('Begin train')
    history = model.fit(
    [X_seq_train, X_numeric_train], y_train,
    validation_split=0.2,
    epochs=20,
    batch_size=32
    )

    loss, accuracy = model.evaluate([X_seq_test, X_numeric_test], y_test)
    print(f"Test Accuracy: {accuracy:.2f}")


