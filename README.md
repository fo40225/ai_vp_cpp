# ai_vp_cpp

```bash
sudo docker build -t aivpcpp:0.0.1 .
sudo docker run -d -p 80:8080 \
--restart unless-stopped \
-e METAMAP_HOST=127.0.0.1 \
-e METAMAP_PORT=5000 \
-e METAMAP_ENDPOINT=/api/api/ \
-e VP_HOST=127.0.0.1 \
-e VP_PORT=80 \
-e VP_ENDPOINT=/VPHandler.ashx \
aivpcpp:0.0.1
```
