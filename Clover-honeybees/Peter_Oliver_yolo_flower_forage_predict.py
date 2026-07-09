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
    Main function to run YOLO flower detection and forage classification on a folder of unlabeled images.
    
    Arguments:
        --flower-path: Path to the trained YOLO flower detection model.
        --forage-path: Path to the trained YOLO forage classification model.
        --unlabeled-dir: Folder with unlabeled images.
        --output-dir: Where to save prediction results.

    Returns:
        None. Saves prediction results in the specified output directory along with a summary CSV file.
    """

    parser = argparse.ArgumentParser(description="YOLO Prediction Script")
    parser.add_argument('--flower-path', type=str, required=True, help='Path to trained YOLO flower detection model')
    parser.add_argument('--forage-path', type=str, required=True, help='Path to trained YOLO forage classification model')
    parser.add_argument('--unlabeled-dir', type=str, required=True, help='Folder with unlabeled images')
    parser.add_argument('--output-dir', type=str, required=True, help='Where to save prediction results')
    args = parser.parse_args()

    FLOWER_PATH = args.flower_path
    FORAGE_PATH = args.forage_path
    UNLABELED_DIR = args.unlabeled_dir
    OUTPUT_DIR = args.output_dir
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load flower detection model
    flower_model = YOLO(FLOWER_PATH)
    flower_model.fuse()
    flower_model.eval()

    # Load forage classification model
    forage_model = YOLO(FORAGE_PATH);
    forage_model.fuse()
    forage_model.eval()
    
    # Move models to device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    flower_model.to(device)
    forage_model.to(device)

    # Get image files
    image_files = [f for f in os.listdir(UNLABELED_DIR) if f.lower().endswith(('.jpg', '.jpeg', '.png'))]

    # Prepare CSV output
    csv_path = os.path.join(OUTPUT_DIR, "summary.csv")
    with open(csv_path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['image', 'datetime', 'num_flowers', 'num_insects'])

        for img_name in image_files:

            # get flower detections
            img_path = os.path.join(UNLABELED_DIR, img_name)
            results = flower_model.predict(img_path, conf=0.75, verbose=False)
            boxes = results[0].boxes.xyxy.cpu().numpy()  # xyxy format
            scores = results[0].boxes.conf.cpu().numpy()
            classes = results[0].boxes.cls.cpu().numpy()
            flower_class_names = results[0].names

            # Load original image once
            orig_img = Image.open(img_path).convert("RGB")
            width, height = orig_img.size

            flower_count = 0
            forage_count = 0

            # Save results to txt file (YOLO format: class x1 y1 x2 y2 flower_conf forage class insect_conf)
            out_path = os.path.join(OUTPUT_DIR, os.path.splitext(img_name)[0] + '.txt')
            with open(out_path, 'w') as f:

                cutouts = []
                coords = []
                for box, score, cls in zip(boxes, scores, classes):
                    if score < 0.75:  # Skip low confidence detections
                        continue
                    flower_count += 1

                    x1, y1, x2, y2 = map(int, box)
                    # Add 20 pixels padding, ensuring we stay within image bounds
                    pad = 20
                    x1_pad = max(x1 - pad, 0)
                    y1_pad = max(y1 - pad, 0)
                    x2_pad = min(x2 + pad, width)
                    y2_pad = min(y2 + pad, height)
                    # Ensure valid crop coordinates
                    if x2_pad <= x1_pad or y2_pad <= y1_pad:
                        continue  # skip invalid crop
                    cutout = orig_img.crop((x1_pad, y1_pad, x2_pad, y2_pad))
                    cutouts.append(cutout)
                    coords.append((x1, y1, x2, y2, score, cls))

                if cutouts:
                    outputs = forage_model.predict(cutouts, verbose=False)
                    predicted_classes = [output.names[output.probs.top1] for output in outputs]
                    confidences = [output.probs.top1conf for output in outputs]

                for i, (x1, y1, x2, y2, score, cls) in enumerate(coords):
                    predicted_class = predicted_classes[i]
                    confidence = confidences[i]
                    flower_class_name = flower_class_names[int(cls)] if flower_class_names and int(cls) < len(flower_class_names) else str(cls)
                    if predicted_class == 'foraging' and confidence >= 0.50:
                        forage_count += 1
                    # Write flower detection and forage classification results to txt file
                    f.write(f"{flower_class_name} {x1} {y1} {x2} {y2} {score:.4f} {predicted_class} {confidence:.4f}\n")

            # Parse datetime from filename
            dt_str = parse_datetime_from_filename(img_name)
            # Write summary row for this image
            csvwriter.writerow([img_name, dt_str, flower_count, forage_count])

    print('All predictions complete!')

if __name__ == "__main__":
    main()
