# ee274-bwt-project

## 0. Getting started
- Create conda environment and install required packages:
    ```
    conda create --name bwt_env python=3.9
    conda activate bwt_env
    ```
- Clone the repo
    ```
    git clone https://github.com/ujinsong23/ee274_bwt_project.git
    cd ee274_bwt_project
    ```
- Install the `scl` package
    ```
    pip install -e . #install the package in a editable mode
    ``` 


`bwt.py`, `compress.py`, `search.py` files in this directory include our implementation of the **efficient BWT** and **searching algorithm** using the FM index.

## 1. `bwt.py` 
- Includes all the necessary functions for our implementation. 
- You may run tests to ensure the the library is installed correctly and all the transforms are working properly.

    ```
    pytest bwt.py
    ``` 

## 2. `compress.py`
- Allows you to compress (from `.txt` to `.bwtz`) or decompress (from `.bwtz` to `.txt`) files.
    ```
    python compress.py [-h] -i INPUT [-d] [-o OUTPUT] [--delimiter DELIMITER]
    ```
- For example,
    - if you want to compress `sherlock_ascii.txt`, 
        ```
        python compress.py -i sherlock_ascii.txt
        ```
    - if you want to decompress `sherlock_ascii.bwtz`, 
        ```
        python compress.py -i sherlock_ascii.bwtz -d
        ```

## 3. `search.py`
- allows you to search for a specific string pattern and count the total occurrences in either the compressed(`.bwtz` format) or the decompressed(`.txt` format) files.
    ```
    python search.py [-h] -i INPUT [-q QUERY] [-d] [-c]
    ```
- For example,
    - if you want to know the occurrences of "Chapter" in **compressed** file `sherlock_ascii.bwtz`, 
        ```
        python search.py -i "sherlock_ascii.bwtz" -q "Chapter" 
        ```
    - if you want to know the occurrences of "Chapter" in decompressed text file `sherlock_ascii.txt`, 
        ```
        python search.py -i "sherlock_ascii.txt" -q "Chapter" -d
        ```
