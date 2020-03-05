#!/bin/bash

python3 -m venv gs2_testing_virtualenv
source gs2_testing_virtualenv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
