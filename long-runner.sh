#!/bin/bash

# Function to generate a random number between 1 and 100
generate_random_number() {
  echo $(( ( RANDOM % 100 )  + 1 ))
}

# Main loop
while true; do
  # Generate a random number
  random_number=$(generate_random_number)
  
  # Output the random number to a file
  echo $random_number >> random_numbers.txt
  
  # Print the generated random number to the console
  echo "Generated random number: $random_number"
  
  # Sleep for 5 seconds
  sleep 10
done
