In order to run this in Docker:

```bash 
heroku container:login
docker build --rm -t chasejohnh/alignment-viewer .
docker run -p 5000:5000 -e PORT=5000 chasejohnh/alignment-viewer:latest
```
