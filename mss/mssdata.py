#import tensorflow as tf 
import os
import pickle


this_dir, this_filename = os.path.split(__file__)

class peakmodel:
	#test_model_path = os.path.join(this_dir, 'model_5_3_14_4_48.h5')
	#test_model = tf.keras.models.load_model(test_model_path)
	Model_file = os.path.join(this_dir, 'rfmodel.pkl')
	rf_model = pickle.load(open(Model_file, 'rb'))

