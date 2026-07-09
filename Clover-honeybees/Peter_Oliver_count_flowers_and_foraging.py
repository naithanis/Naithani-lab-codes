

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
This script processes a directory of .txt files containing flower detections and
foraging classifications. It aggregates the counts of flowers and unique foraging 
interactions over 10-minute intervals and outputs a summary CSV file.
"""

import os
import sys
import glob
import csv
from collections import Counter, deque
from datetime import datetime, timedelta

def parse_line(line):
    # Parse a line into its components
    parts = line.strip().split()
    if len(parts) != 8:
        return None
    flower_class = parts[0]
    bbox = list(map(int, parts[1:5]))
    flower_conf = float(parts[5])
    foraging_pred = parts[6]
    foraging_conf = float(parts[7])
    return {
        'flower_class': flower_class,
        'bbox': bbox,
        'flower_conf': flower_conf,
        'foraging_pred': foraging_pred,
        'foraging_conf': foraging_conf
    }

def process_file(filepath, flower_conf_threshold=0.9,foraging_conf_threshold=0.98):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    flower_count = 0
    foraging_count = 0
    high_conf_bboxes = []
    for line in lines:
        parsed = parse_line(line)
        if not parsed:
            continue
        flower_count += 1
        if parsed['flower_conf'] > flower_conf_threshold and parsed['foraging_pred'] == 'foraging' and parsed['foraging_conf'] > foraging_conf_threshold:
            foraging_count += 1
            high_conf_bboxes.append(parsed['bbox'])
    return flower_count, foraging_count, high_conf_bboxes

def parse_filename(filename):
    # Expects format: Camera_YearMonthDay_HourMinuteSecond.txt
    base = os.path.basename(filename)
    parts = base.split('_')
    if len(parts) < 3:
        return None
    date_str = parts[1]
    time_str = parts[2].split('.')[0]
    dt = datetime.strptime(date_str + time_str, "%Y%m%d%H%M%S")
    return dt

def bbox_distance(b1, b2):
    # Simple sum of absolute differences for all coordinates
    return sum(abs(a - b) for a, b in zip(b1, b2))

def get_interval_start(dt):
    # Find the nearest previous 5, 15, 25, 35, 45, or 55 minute mark
    minute = dt.minute
    interval_minute = max([m for m in [5, 15, 25, 35, 45, 55] if m <= minute], default=55)
    if interval_minute == 55 and minute < 5:
        dt = dt - timedelta(hours=1)
    return dt.replace(minute=interval_minute, second=0, microsecond=0)

def interval_midpoint(start):
    return start + timedelta(minutes=5)

def main():
    """
    Main function to process a directory of .txt files containing flower and foraging predictions.
    Aggregates counts over 10-minute intervals and writes a summary CSV file.

    Arguments:
        directory: str, path to the directory containing .txt files for each image prediction.
        n_prev: int, number of previous image results to consider for unique foraging counts.

    Returns:
        None. Writes a summary CSV file in the specified directory.
    """

    if len(sys.argv) < 2:
        print("Usage: python count_flowers_and_foraging.py <directory> [n_prev]")
        sys.exit(1)
    directory = sys.argv[1]
    n_prev = int(sys.argv[2]) if len(sys.argv) > 2 else 2  # Now configurable via CLI
    pixel_threshold = 80  # Max sum(abs(coord diff)) to consider same insect (20 pixels on all sides)

     # Extract MonthDayYear and CameraNumber from directory name
    dir_base = os.path.basename(os.path.normpath(directory))
    dir_parts = dir_base.split('_')
    if len(dir_parts) >= 3:
        month_day_year = dir_parts[1]
        camera_number = dir_parts[2]
        out_csv = os.path.join(directory, f"{camera_number}_{month_day_year}_counts.csv")
    else:
        out_csv = os.path.join(directory, "interval_summary.csv")

    txt_files = sorted(glob.glob(os.path.join(directory, "*.txt")))
    if not txt_files:
        print("No .txt files found.")
        sys.exit(1)

     # Determine first interval start
    first_dt = parse_filename(txt_files[0])
    if not first_dt:
        print("Could not parse datetime from first file.")
        sys.exit(1)
    current_interval_start = get_interval_start(first_dt)
    current_interval_end = current_interval_start + timedelta(minutes=10)

    interval_flower_counts = []
    interval_unique_foraging_count = 0
    recent_foraging_bboxes = deque(maxlen=n_prev)  # Store bboxes from last n_prev intervals

    with open(out_csv, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["interval_midpoint", "mode_flower_count", "unique_foraging_count"])

        for txt_file in txt_files:
            dt = parse_filename(txt_file)
            if not dt:
                continue

            # If file is outside current interval, aggregate and write, then reset
            while dt >= current_interval_end:
                if interval_flower_counts:
                    mode_bbox_count = Counter(interval_flower_counts).most_common(1)[0][0]
                    midpoint = interval_midpoint(current_interval_start)
                    writer.writerow([midpoint.strftime("%Y-%m-%d %H:%M:%S"), mode_bbox_count, interval_unique_foraging_count])
                # Prepare for next interval
                current_interval_start += timedelta(minutes=10)
                current_interval_end = current_interval_start + timedelta(minutes=10)
                interval_flower_counts = []
                interval_unique_foraging_count = 0

            flower_count, foraging_count, high_conf_bboxes = process_file(txt_file)
            interval_flower_counts.append(flower_count)
            # Check for unique foraging interactions
            for bbox in high_conf_bboxes:
                new_bbox = True
                for prev_foraging in recent_foraging_bboxes:
                    if any(bbox_distance(bbox, prev_bbox) <= pixel_threshold for prev_bbox in prev_foraging):
                        new_bbox = False
                        break
                if new_bbox:
                    interval_unique_foraging_count += 1
            recent_foraging_bboxes.append(high_conf_bboxes)
            if len(recent_foraging_bboxes) > n_prev:
                recent_foraging_bboxes.popleft()

        # Write last interval if any files remain
        if interval_flower_counts:
            mode_bbox_count = Counter(interval_flower_counts).most_common(1)[0][0]
            midpoint = interval_midpoint(current_interval_start)
            writer.writerow([midpoint.strftime("%Y-%m-%d %H:%M:%S"), mode_bbox_count, interval_unique_foraging_count])

    print(f"Interval summary written to {out_csv}")

if __name__ == "__main__":
    main()
