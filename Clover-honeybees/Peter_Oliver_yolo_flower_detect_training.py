
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
Train a YOLO model for flower detection using transfer learning.

This script uses the Ultralytics YOLO library to train a YOLO model for 
clover flower detection. It assumes that you have a labeled dataset of 
cloverflower images and corresponding annotations in YOLO format. Relative 
paths to the dataset are expected to be structured as follows:

./data/labeled/
    images/
    labels/
    classes.txt

The script automatically splits the dataset into training and validation sets.
The final trained weights are saved in a ./models directory.
"""

import os
from ultralytics import YOLO
from ultralytics.data.split import autosplit

# Paths
DATA_DIR = './data/labeled'
IMAGES_DIR = os.path.join(DATA_DIR, 'images')
LABELS_DIR = os.path.join(DATA_DIR, 'labels')
CLASSES_FILE = os.path.join(DATA_DIR, 'classes.txt')

# Autosplit: creates train/val split text files in LABELS_DIR
print('Splitting dataset with autosplit...')
autosplit(IMAGES_DIR, weights=(0.8, 0.2, 0.0))  # 80% train, 20% val, 0% test

# Prepare YOLO dataset config
data_yaml = f"""
path: {DATA_DIR}
train: autosplit_train.txt
val: autosplit_val.txt
test: autosplit_test.txt
nc: {len(open(CLASSES_FILE).readlines())}
names:
"""
with open(CLASSES_FILE, 'r') as f:
    for idx, name in enumerate(f):
        data_yaml += f"  {idx}: {name.strip()}\n"

with open('yolo_flower_detect_training.yaml', 'w') as f:
    f.write(data_yaml)

# Load pretrained YOLO model
model = YOLO('yolo11l.pt')

# Train with transfer learning
model.train(
    data='yolo_flower_detect_training.yaml',
    imgsz=2048,
    epochs=50,
    batch=4,
    workers=8,
    device='cuda',
    name='yolo11_flower_detect',
    conf=0.5
)

print('Training complete!')

# Save the trained model to the models directory
MODELS_DIR = './models'
os.makedirs(MODELS_DIR, exist_ok=True)
model_path = os.path.join(MODELS_DIR, 'yolo11_flower_detect.pt')
model.save(model_path)
print(f'Model saved to {model_path}')
