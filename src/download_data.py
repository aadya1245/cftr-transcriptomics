"""Download missing public inputs and verify pinned SHA256 checksums."""
import hashlib,json,urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    for item in json.loads((ROOT/'data/manifest.json').read_text()):
        path=ROOT/item['path']
        if not path.exists():
            payload=urllib.request.urlopen(item['url'],timeout=120).read()
            if hashlib.sha256(payload).hexdigest()!=item['sha256']:
                raise ValueError('Downloaded source has changed: '+item['path'])
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes(payload)
        if hashlib.sha256(path.read_bytes()).hexdigest()!=item['sha256']:
            raise ValueError('Checksum mismatch: '+item['path'])
        print('Verified',item['path'])
if __name__=='__main__': main()
