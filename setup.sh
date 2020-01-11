#!/usr/bin/env bash

echo "Setting up for HW1..."

if ! hash wget 2>/dev/null; then
  echo "Installing wget..."
  sudo apt-get install wget
fi

echo "Installing python requirements (including gdown)..."
python3 -m pip install -r requirements.txt

echo "Fetching PDF instructions..."
wget https://faculty.washington.edu/kutz/582hw1.pdf
echo "Fetching test data..."
gdown https://drive.google.com/uc?id=1mQnKu6GO8x130tiaoj0dDty97PgtGzpF 
