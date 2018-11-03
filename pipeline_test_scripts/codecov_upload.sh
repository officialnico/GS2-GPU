#!/bin/bash

if [ $? -eq 0 ]; then bash <(curl -s https://codecov.io/bash) -t "$CODECOV_TOKEN"; fi

