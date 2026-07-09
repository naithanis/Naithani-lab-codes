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


import os
import argparse
from ultralytics import YOLO
from PIL import Image
import torch
import torchvision
from torchvision import transforms
import csv
from datetime import datetime

def parse_datetime_from_filename(filename):
    """
    Parse datetime from filename in the format CameraNumber_YEARMONTHDAY_HOURMINUTESECOND.jpg
    Example: C1_20240607_153012.jpg
    """
    base = os.path.splitext(os.path.basename(filename))[0]
    parts = base.split('_')
    if len(parts) >= 3:
        date_str = parts[1]
        time_str = parts[2]
        try:
            dt = datetime.strptime(date_str + time_str, "%Y%m%d%H%M%S")
            return dt.strftime("%Y-%m-%d %H:%M:%S")
        except Exception:
            return ""
    return ""

def main():
    """
    Main function to run YOLO flower detection and insect classification on a folder of unlabeled images.
    
    Arguments:
        --flower-path: Path to the trained YOLO flower detection model.
        --insect-path: Path to the trained regnet insect classification model.
                       (from https://github.com/ShivaniChiranjeevi/Insect-Classifier/tree/main)
        --unlabeled-dir: Folder with unlabeled images.
        --output-dir: Where to save prediction results.

    Returns:
        None. Saves prediction results in the specified output directory with the following subdirectories:
            --OUTPUT_BEE: Contains cutouts of flowers classified as having a foraging insect on them
            --OUTPUT_NOBEE: Contains cutouts of flowers classified as not having a foraging insect on them
    """

    parser = argparse.ArgumentParser(description="YOLO Prediction Script")
    parser.add_argument('--flower-path', type=str, required=True, help='Path to trained YOLO flower detection model')
    parser.add_argument('--insect-path', type=str, required=True, help='Path to trained regnet insect classification model')
    parser.add_argument('--unlabeled-dir', type=str, required=True, help='Folder with unlabeled images')
    parser.add_argument('--output-dir', type=str, required=True, help='Where to save prediction results')
    args = parser.parse_args()

    FLOWER_PATH = args.flower_path
    INSECT_PATH = args.insect_path
    UNLABELED_DIR = args.unlabeled_dir
    OUTPUT_DIR = args.output_dir
    OUTPUT_BEE = os.path.join(OUTPUT_DIR, "foraging")
    OUTPUT_NOBEE = os.path.join(OUTPUT_DIR, "noforaging")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(OUTPUT_BEE, exist_ok=True)
    os.makedirs(OUTPUT_NOBEE, exist_ok=True)

    # Load flower detection model
    flower_model = YOLO(FLOWER_PATH)
    flower_model.fuse()
    flower_model.eval()
    
    # Load insect classification model
    insect_model = torchvision.models.regnet_y_32gf()
    # Set the correct output layer before loading weights
    insect_model.fc = torch.nn.Linear(insect_model.fc.in_features, 2526)
    
    # load insect model weights
    checkpoint = torch.load(INSECT_PATH, map_location=torch.device('cpu'), weights_only=False)
    insect_model.load_state_dict(checkpoint['model'], strict=True)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    insect_model.eval()
    
    # Move models to device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    flower_model.to(device)
    insect_model.to(device)

    # image preprocessing for cutouts
    preprocess = transforms.Compose([
        transforms.Resize(256),
        transforms.CenterCrop(224),
        transforms.ToTensor(),
        transforms.ConvertImageDtype(torch.float),
        transforms.Normalize(
            mean=[0.485, 0.456, 0.406],
            std=[0.229, 0.224, 0.225]
        ),
    ])

    # Get image files
    image_files = [f for f in os.listdir(UNLABELED_DIR) if f.lower().endswith(('.jpg', '.jpeg', '.png'))]

    for img_name in image_files:
        img_path = os.path.join(UNLABELED_DIR, img_name)
        results = flower_model.predict(img_path, conf=0.90)
        boxes = results[0].boxes.xyxy.cpu().numpy()  # xyxy format
        scores = results[0].boxes.conf.cpu().numpy()
        classes = results[0].boxes.cls.cpu().numpy()

        # Load original image once
        orig_img = Image.open(img_path).convert("RGB")
        width, height = orig_img.size

        cutouts = []
        coords = []
        for box, score, cls in zip(boxes, scores, classes):
            if score < 0.90:
                continue
            x1, y1, x2, y2 = map(int, box)
            pad = 20
            x1_pad = max(x1 - pad, 0)
            y1_pad = max(y1 - pad, 0)
            x2_pad = min(x2 + pad, width)
            y2_pad = min(y2 + pad, height)
            if x2_pad <= x1_pad or y2_pad <= y1_pad:
                continue  # skip invalid crop
            cutout = orig_img.crop((x1_pad, y1_pad, x2_pad, y2_pad))
            cutouts.append(cutout)
            coords.append((img_name, i if 'i' in locals() else 0, x1, y1, x2, y2, score, cls))

        if cutouts:
            input_batch = torch.stack([preprocess(c) for c in cutouts]).to(device)
            with torch.no_grad():
                outputs = insect_model(input_batch)
                probabilities = torch.nn.functional.softmax(outputs, dim=1)
                predicted_classes = torch.argmax(probabilities, dim=1).cpu().numpy()
                confidences = probabilities.max(dim=1).values.cpu().numpy()

        # Count files in OUTPUT_BEE and OUTPUT_NOBEE (all files are .jpg)
        bee_files = len(os.listdir(OUTPUT_BEE))
        nobee_files = len(os.listdir(OUTPUT_NOBEE))

        for i, (img_name, _, x1, y1, x2, y2, score, cls) in enumerate(coords):
            predicted_class = predicted_classes[i]
            confidence = confidences[i]
            cutout_filename = f"{os.path.splitext(img_name)[0]}_flower{i}.jpg"
            if confidence > 0.19 and predicted_class in [498,501,504]: # classes matching honeybee and related species
                cutout_path = os.path.join(OUTPUT_BEE, cutout_filename)
                cutouts[i].save(cutout_path)
                bee_files += 1
            elif confidence < 0.1 and nobee_files < bee_files:
                cutout_path = os.path.join(OUTPUT_NOBEE, cutout_filename)
                cutouts[i].save(cutout_path)
                nobee_files += 1

    print('All predictions complete!')

if __name__ == "__main__":
    main()
