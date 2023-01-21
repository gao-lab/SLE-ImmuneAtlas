# reappearance  R env in Rstudio

### Push the docker

```bash
docker push rocker/rstudio:4.0.3
docker run -d --name sle_rstudio --hostname rstudio -P -p 8889:8787 -v /rd2/user/xiacr:/data rocker/rstudio:4.0.3
```
NOTE: if you face problem when login, please clean the cache in browser


### Install env correlated packages
Run in Rstudio
```R
install.packages('renv')
install.packages('pak')
```

### Install main packages

```R
setwd('/data/sle/env')
renv::
```s
