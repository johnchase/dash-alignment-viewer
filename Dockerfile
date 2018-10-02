FROM heroku/miniconda:3

# Grab requirements.txt.
ADD ./webapp/requirements.txt /tmp/requirements.txt

# Install dependencies
RUN pip install -r /tmp/requirements.txt

# Add our code
ADD ./webapp /opt/webapp/
WORKDIR /opt/webapp


RUN apt-get update && \
    apt-get install libgl1-mesa-glx -y

RUN conda install scikit-bio
RUN useradd -m myuser
USER myuser

CMD gunicorn --bind 0.0.0.0:$PORT wsgi
