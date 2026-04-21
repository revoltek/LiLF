#!/usr/bin/env python3
# Generic file download with retry and check for length

import os
from time import sleep
import requests
from requests.packages.urllib3.exceptions import InsecureRequestWarning

requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

def download_file(url, filename, login=None, password=None,macaroon=None):

    downloaded=False
    while not downloaded:
        try:
            macaroon=list(macaroon.values())[0]
        except:
            # In this case, the macaroon was already stored as a string, so nothing has to change
            pass
        connected=False
        tries=1
        while not connected:
            try:
                #print('Opening connection')
                if (login is not None) and (password is not None):
                    response = requests.get(url+f'?authz={macaroon}', stream=True, verify=False, timeout=60,
                            auth=(login, password),headers={'Authorization':f'bearer {macaroon}'})
                else:
                    response = requests.get(url, stream=True, verify=True, timeout=60)
                    print('nooo')

                if response.status_code!=200:
                    raise RuntimeError('Code was %i' % response.status_code)
                if 'Content-Length' in response.headers.keys():
                    esize = int(response.headers['Content-Length'])
                else: esize = None
            except requests.exceptions.ConnectionError:
                print('Downloader -- Connection error! sleeping 30 seconds before retry (%i/10)...' % tries)
                sleep(30)
            except (requests.exceptions.Timeout,requests.exceptions.ReadTimeout):
                print('Downloader -- Timeout! sleeping 30 seconds before retry (%i/10)...' % tries)
                sleep(30)
            except RuntimeError:
                print('Downloader -- Download error (%i/10)...' % tries)
                sleep(30)
            else:
                connected=True
            if tries > 10:
                print('Downloader -- Giving up downlaoding the file...')
                return False # false
            tries += 1
        try:
            #print('Downloading %i bytes' % esize)
            with open(filename, 'wb') as fd:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        fd.write(chunk)
            fsize=os.path.getsize(filename)
            if esize!=fsize and esize is not None:
                print('Downloader -- Download incomplete (expected %i, got %i)! Retrying...' % (esize, fsize))
            elif esize is not None:
                print('Downloader -- Download successful, all %i bytes received' % (fsize))
                downloaded=True
            else:
                print('Downloader -- Download successful, %i bytes (unknown size)' % (fsize))
                downloaded=True

        except (requests.exceptions.ConnectionError,requests.exceptions.Timeout,requests.exceptions.ChunkedEncodingError):
            print('Downloader -- Connection error! sleeping 30 seconds before retry...')
            sleep(30) # back to the connection

    del response
    return downloaded

if __name__=='__main__':
    import sys
    download_file(sys.argv[1],sys.argv[2])
