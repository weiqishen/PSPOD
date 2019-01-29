import os

def Fatal_Error(err_str, err_code=1):
    print(err_str)
    os._exit(err_code)