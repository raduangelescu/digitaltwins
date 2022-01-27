import py7zr
import requests
import os

def download_file_from_google_drive(id, destination):
    URL = "https://docs.google.com/uc?export=download"

    session = requests.Session()

    response = session.get(URL, params = { 'id' : id }, stream = True)
    token = get_confirm_token(response)

    if token:
        params = { 'id' : id, 'confirm' : token }
        response = session.get(URL, params = params, stream = True)

    save_response_content(response, destination)    

def get_confirm_token(response):
    for key, value in response.cookies.items():
        if key.startswith('download_warning'):
            return value

    return None

def save_response_content(response, destination):
    CHUNK_SIZE = 32768

    with open(destination, "wb") as f:
        for chunk in response.iter_content(CHUNK_SIZE):
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)
file_id = "1qIR-tpIIGM8LRit2A5BJZdxRLHUyrq9q"
print("downloading data file")
tmp_download_filename = 'data_tmp.7z'
download_file_from_google_drive(file_id, tmp_download_filename)
print('finished downloading, unpacking data')
with py7zr.SevenZipFile(tmp_download_filename, mode='r') as z:
    z.extractall('data')
print('finished unpacking data, cleaning temp files')
os.remove(tmp_download_filename)
print('done')