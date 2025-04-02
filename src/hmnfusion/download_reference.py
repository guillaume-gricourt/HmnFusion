import hashlib
import logging
import urllib.parse
from typing import Dict

import requests


class DownloadZenodo(object):
    ZENODO_API = "https://zenodo.org/api/"

    def __init__(self, id: str, params: Dict[str, str]) -> None:
        self.id = id
        self.params = params

    def show_files(self) -> Dict[str, str]:
        url = urllib.parse.urljoin(
            self.ZENODO_API, "deposit/depositions/%s/files" % (self.id,)
        )
        logging.info("url: %s" % (url,))
        r = requests.get(url, params=self.params)
        self.check_request(request=r)
        return r.json()

    def download_file(self, url: str, path: str) -> None:
        try:
            response = requests.get(url, stream=True, params=self.params)
            self.check_request(request=response)
            with open(path, "wb") as fod:
                for chunk in response.iter_content(chunk_size=8192):
                    fod.write(chunk)
            response.close()
        except Exception as e:
            raise ValueError(str(e))

    @classmethod
    def check_request(cls, request) -> None:
        if request.status_code > 202:
            raise ValueError(request.text)

    @classmethod
    def check_checksum(cls, checksum: str, path: str) -> None:
        hash_md5 = hashlib.md5()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        if checksum != hash_md5.hexdigest():
            raise ValueError("Checksum computed is not conform for file: %s" % (path,))

    @classmethod
    def load_token(cls, path: str) -> str:
        access_token = ""
        with open(path) as fid:
            data = fid.read()
            access_token = data.replace("\n", "")
            assert access_token is str
            assert len(access_token) > 1
        return access_token
