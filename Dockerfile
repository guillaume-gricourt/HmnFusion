FROM python:3.9

COPY *whl /opt/
RUN pip install --no-cache-dir /opt/*whl


wget https://github.com/OpenGene/GeneFuse/archive/refs/tags/v0.8.0.tar.gz -O GeneFuse.tar.gz

# build
cd genefuse
make
sudo make install

git clone --recursive https://github.com/arq5x/lumpy-sv.git
wget https://github.com/arq5x/lumpy-sv/archive/refs/tags/v0.3.0.tar.gz -O Lumpy.tar.gz
cd lumpy-sv
make
cp bin/* /usr/local/bin/.

ENTRYPOINT ["hmnfusion"]
