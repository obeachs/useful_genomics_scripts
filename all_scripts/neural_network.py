import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.model_selection import train_test_split
from tensorflow.keras.layers import Conv1D, Dense, MaxPooling1D, Flatten
from tensorflow.keras.models import Sequential
from sklearn.metrics import confusion_matrix
import itertools
import tensorflow.keras.backend as K
SEQUENCES_URL = 'https://raw.githubusercontent.com/abidlabs/deep-learning-genomics-primer/master/sequences.txt'
sequences = requests.get(SEQUENCES_URL).text.split('\n')
sequences = list(filter(None,sequences))
print(sequences)

df = pd.DataFrame(sequences, index=np.arange(1,len(sequences)+1),columns=['Sequences'])
#print(df.head())

integer_encoder = LabelEncoder()
one_hot_encoder = OneHotEncoder(categories = 'auto')
input_features = []
count = 0
#print(integer_encoder.fit_transform(list(sequences[1])))
for sequence in sequences:
    '''This is turning each sequence in the sequence list into numbers (0-3) for
     0=A,C=1,G=2,T=3'''
    integer_encoded = integer_encoder.fit_transform(list(sequence))
    '''Now using numpy to create an array, one column, one row'''
    integer_encoded = np.array(integer_encoded.reshape(-1,1))
    '''Another possibility to convert categorical features to features that can be used with scikit-learn estimators is to use a one-of-K,
     also known as one-hot or dummy encoding. This type
     of encoding can be obtained with the OneHotEncoder, which transforms each
    categorical feature with n_categories possible values into n_categories binary features, with one of them 1, and all others 0.'''
    one_hot_encoded = one_hot_encoder.fit_transform(integer_encoded)
    input_features.append(one_hot_encoded.toarray())
shape = input_features[1].shape
for array in input_features:
    if array.shape != shape:
        print(array.shape)
print(shape)
print(input_features)
#print(count)
np.set_printoptions(threshold = 40)
input_features = np.stack(input_features)
print('Example sequences\n----------------')
print('DNA sequence 1:\n',sequences[0][:10],'....',sequences[0][-10:])
'''.T in this case is just transposing the matrix'''
print('One hot encoding of Sequences 1:\n',input_features[0].T)


LABELS_URL = 'https://raw.githubusercontent.com/abidlabs/deep-learning-genomics-primer/master/labels.txt'

labels = requests.get(LABELS_URL).text.split('\n')
labels = list(filter(None, labels))  # removes empty sequences

one_hot_encoder = OneHotEncoder(categories='auto')
labels = np.array(labels).reshape(-1, 1)
input_labels = one_hot_encoder.fit_transform(labels).toarray()

print('Labels:\n',labels.T)
print('One-hot encoded labels:\n',input_labels.T)

'''This part is taking advantage of the scikit-learn module, it is splitting the
input features and labels that were made using one_hot_encoder into random subsets
made up of test and training inputs '''
train_features, test_features, train_labels, test_labels = train_test_split(
    input_features, input_labels, test_size=0.25, random_state=42)

'''Now we are using the tensorflow keras modules to create a model for our CNN
transcription factor binding site predictor'''
print(train_features)
model = Sequential()
'''32 filters with base size of 12'''
model.add(Conv1D(filters=32, kernel_size=12,
                 input_shape=(train_features.shape[1], 4)))
model.add(MaxPooling1D(pool_size=4))
model.add(Flatten())
model.add(Dense(16, activation='relu'))
model.add(Dense(2, activation='softmax'))

model.compile(loss='binary_crossentropy', optimizer='adam',
              metrics=['binary_accuracy'])
model.summary()

'''model.fit() is the method we run to run the neural network. epochs variable
defines how many iterations it does - we want to make sure that the network runs
just enough that it gets a good prediction rate but does not dip too far below its
peak'''
history = model.fit(train_features, train_labels,
                    epochs=50, verbose=0, validation_split=0.25)

#plt.figure()
#plt.plot(history.history['loss'])
#plt.plot(history.history['val_loss'])
#plt.title('model loss')
#plt.ylabel('loss')
#plt.xlabel('epoch')
#plt.legend(['train', 'validation'])
#plt.show()

#plt.figure()
#plt.plot(history.history['binary_accuracy'])
#plt.plot(history.history['val_binary_accuracy'])
#plt.title('model accuracy')
#plt.ylabel('accuracy')
#plt.xlabel('epoch')
#plt.legend(['train', 'validation'])
#plt.show()

predicted_labels = model.predict(np.stack(test_features))
cm = confusion_matrix(np.argmax(test_labels, axis=1),
                      np.argmax(predicted_labels, axis=1))
print('Confusion matrix:\n',cm)

cm = cm.astype('float') / cm.sum(axis = 1)[:, np.newaxis]

#plt.imshow(cm, cmap=plt.cm.Blues)
#plt.title('Normalized confusion matrix')
#plt.colorbar()
#plt.xlabel('True label')
#plt.ylabel('Predicted label')
#plt.xticks([0, 1]); plt.yticks([0, 1])
#plt.grid('off')
for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
    plt.text(j, i, format(cm[i, j], '.2f'),
             horizontalalignment='center',
             color='white' if cm[i, j] > 0.5 else 'black')
plt.show()

def compute_salient_bases(model, x):
  input_tensors = [model.input]
  gradients = model.optimizer.get_gradients(model.output[0][1], model.input)
  compute_gradients = K.function(inputs = input_tensors, outputs = gradients)

  x_value = np.expand_dims(x, axis=0)
  gradients = compute_gradients([x_value])[0][0]
  sal = np.clip(np.sum(np.multiply(gradients,x), axis=1),a_min=0, a_max=None)
  return sal


sequence_index = 1999  # You can change this to compute the gradient for a different example. But if so, change the coloring below as well.
sal = compute_salient_bases(model, input_features[sequence_index])

plt.figure(figsize=[16,5])
barlist = plt.bar(np.arange(len(sal)), sal)
[barlist[i].set_color('C2') for i in range(5,17)]  # Change the coloring here if you change the sequence index.
plt.xlabel('Bases')
plt.ylabel('Magnitude of saliency values')
plt.xticks(np.arange(len(sal)), list(sequences[sequence_index]));
plt.title('Saliency map for bases in one of the positive sequences'
          ' (green indicates the actual bases in motif)');
plt.show()
