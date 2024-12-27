import pandas as pd
import numpy as np
from pyfaidx import Fasta
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Conv1D, GlobalMaxPooling1D, Dense, Concatenate, Dropout
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

# Enable mixed precision for performance
tf.keras.mixed_precision.set_global_policy('mixed_float16')

# Helper function for one-hot encoding
def prepare_motif_encodings(sequence, max_len):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    encoded = np.zeros((max_len, 4))  # Fixed max_len for uniformity
    for i, base in enumerate(sequence[:max_len]):
        if base in mapping:
            encoded[i, mapping[base]] = 1
    return encoded

if __name__ == "__main__":
    # Load the dataset in chunks to reduce memory usage
    chunks = pd.read_csv('data/expression_tr_annotations.txt', sep='\t', chunksize=100000)

    processed_chunks = []
    fasta = Fasta('data/hg38.ml.fa')
    numeric_features = []  # Placeholder: Define based on your columns

    # Process each chunk
    for chunk in chunks:
        # Prepare sequences and labels
        chunk['sequence'] = chunk.apply(
            lambda row: str(fasta[row['chr']][row['ann_start'] - 1 : row['ann_end']]), axis=1
        )
        chunk['disease_binary'] = chunk['disease'].apply(lambda x: 1 if isinstance(x, str) else 0)
        processed_chunks.append(chunk)

    # Combine processed chunks
    df = pd.concat(processed_chunks)

    # One-hot encode sequences
    max_len = (df['ann_end'] - df['ann_start']).max()
    sequence_encodings = np.array([prepare_motif_encodings(seq, max_len) for seq in df['sequence']])

    breakpoint()
    # Normalize numeric features
    scaler = StandardScaler()
    df[numeric_features] = scaler.fit_transform(df[numeric_features])

    # Split data
    X_seq = sequence_encodings
    X_numeric = df[numeric_features].values
    y = np.array(df['disease_binary'])
    X_seq_train, X_seq_test, X_numeric_train, X_numeric_test, y_train, y_test = train_test_split(
        X_seq, X_numeric, y, test_size=0.2, random_state=42
    )

    # Define model
    seq_input = Input(shape=(max_len, 4), name="sequence_input")
    x = Conv1D(filters=32, kernel_size=10, activation="relu")(seq_input)
    x = GlobalMaxPooling1D()(x)

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

    # Train model
    history = model.fit(
        [X_seq_train, X_numeric_train], y_train,
        validation_split=0.2,
        epochs=20,
        batch_size=32
    )

    # Evaluate model
    loss, accuracy = model.evaluate([X_seq_test, X_numeric_test], y_test)
    print(f"Test Accuracy: {accuracy:.2f}")
