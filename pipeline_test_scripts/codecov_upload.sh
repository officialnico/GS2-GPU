#!/bin/bash

if [ $? -eq 0 ]; then bash <(curl -s https://codecov.io/bash) -t "a848e445-4d00-4aa9-8d00-e9950fcf29f2"; fi

