# microstructure-homogenization

## Installation Instructions

*************************************These instructions assume Python has been installed.************************************* 

1. Download the .zip file from https://github.com/revvu/microstructure-homogenization.git using the “CODE” dropdown. 
2. Move the file to the desired location and unzip the file.
3. Open Windows PowerShell (or Terminal if on another OS). Change directory to the `microstructure-homogenization-main` folder using “cd [location]”
4. Create a new python virtual environment using the command below:
    
    ```powershell
    python -m venv microstructure-homogenization-env
    ```
    
5. Activate the virtual environment using the command below:
    
    For Windows:
    
    ```powershell
    microstructure-homogenization-env\Scripts\activate.bat
    ```
    
    Unix or MacOS:
    
    ```powershell
    source tutorial-env/bin/activate
    ```
    
6. Make sure the folder you are in contains requirements.txt. Install all dependencies using:
    
    ```powershell
    pip install -r requirements.txt
    ```
    
7. Run the application using:
    
    ```powershell
    python flaskmicrostructures.py
    ```
    

The application will be visible at http://127.0.0.1:5000/.
