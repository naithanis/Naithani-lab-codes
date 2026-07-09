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

import splitfolders

# Define the input and output directories
input_folder = './data/classified/cutouts'
output_folder = './data/classified/split'

# Perform the split (80% train, 20% val, 0% test)
splitfolders.ratio(
    input_folder,
    output=output_folder,
    seed=42,
    ratio=(0.8, 0.2, 0.0),
    group_prefix=None,
    move=False
)
