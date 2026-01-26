#!/bin/bash

# Run pytest with coverage
python -m pytest pytest/unit/ -v --cov=bioutils_collection --cov-report=html --cov-report=term
