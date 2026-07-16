# Honeybee Clover Foraging Identification

This project is designed to identify honeybee foraging activity on clover flowers using a combination of YOLO object detection and image classification models. Flower counts and foraging events are recorded over time, given a time series of camera images.

## Requirements

This project depends on Python 3.12, CUDA 12.8, and the libraries listed in `requirements.txt`. You can install the required libraries using pip:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

The `yolo_flower_forage_predict.py` script uses two Ultralytics YOLOv11 models to first detect clover flowers in images and then classify the detected flowers as being actively foraged or not. The script takes in arguments for the flower detection model path, the foraging classification model path, the input image directory, and the output directory for saving results:

```bash
python yolo_flower_forage_predict.py --flower_model_path <path_to_flower_model> --forage_model_path <path_to_forage_model> --input_dir <path_to_input_images> --output_dir <path_to_output_results>
```

The output is a CSV file containing bounding box coordinates for detected flowers, their classification as foraged or not, and the corresponding image filenames.

You can download pre-trained model weights [here](https://oregonstate.box.com/s/11uck8e1102ymlw2mfw0gqqbzb5xgz97). Alternatively, you can train your own models using `yolo_flower_detect_train.py` for flower detection and `yolo_forage_classify_train.py` for foraging classification. For more details on training YOLO models, refer to the [Ultralytics YOLOv11 documentation](https://docs.ultralytics.com/models/yolo11/).

For initial insect identification, this project used a pre-trained RegNet model from the InsectNet project ([Insect-Classifier](https://github.com/ShivaniChiranjeevi/Insect-Classifier/tree/main) on GitHub).

## References

[1] The Ultralytics YOLOv11 components of this repository are adapted from Ultralytics documentation: [https://docs.ultralytics.com/models/yolo11/](https://docs.ultralytics.com/models/yolo11/)

[2] Shivani Chiranjeevi, Mojdeh Saadati, Zi K Deng, Jayanth Koushik, Talukder Z Jubery, Daren S Mueller, Matthew O’Neal, Nirav Merchant, Aarti Singh, Asheesh K Singh, Soumik Sarkar, Arti Singh, Baskar Ganapathysubramanian. (2025) [InsectNet: Real-time identification of insects using an end-to-end machine learning pipeline](https://doi.org/10.1093/pnasnexus/pgae575). PNAS Nexus, Volume 4, Issue 1, pgae575.

