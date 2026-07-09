# Copyright 2026 Peter Oliver, Yue Zhang and Sushma Naithani
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
Train a YOLO model for foraging classification using transfer learning.

This script uses the Ultralytics YOLO library to train a YOLO model for
insect foraging classification. It assumes that you have a labeled dataset
of insect images and organized by class. Relative paths to the dataset are 
expected to be structured with a train/val/test split as follows:

./data/classified/split/
	train/
		foraging/
        noforaging/
	val/
		foraging/
		noforaging/
	test/
		foraging/
        noforaging/
    
The final trained weights are saved in a ./models directory.
"""

import os
from ultralytics import YOLO

# Paths
DATA_DIR = './data/classified/split'
PRETRAINED_MODEL = 'yolo11l-cls.pt'
EPOCHS = 50
BATCH = 32
IMG_SIZE = 224

# Load model (YOLOv8/YOLOv11 for classification)
model = YOLO(PRETRAINED_MODEL)

# Train with transfer learning
model.train(
	data=DATA_DIR,
	epochs=EPOCHS,
	imgsz=IMG_SIZE,
	batch=BATCH,
    workers=8,
    device='cuda',
	name='yolo11_forage_classify'
)

print('Training complete!')

# Save the trained model to the models directory
MODELS_DIR = './models'
os.makedirs(MODELS_DIR, exist_ok=True)
model_path = os.path.join(MODELS_DIR, 'yolo11_forage_classify.pt')
model.save(model_path)
print(f'Model saved to {model_path}')
