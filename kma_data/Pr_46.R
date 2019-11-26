setwd("C:\\Users\\UOS\\Documents\\GITHUB\\gev\\kma_data")
rm( list = ls()); gc(reset = T)
load("kma_data.Rdata")

if(!require(dplyr)){ install.packages('dplyr')}; require(dplyr)
if(!require(data.table)){ install.packages('data.table')}; require(data.table)
if(!require(maps)){ install.packages('maps')}; require(maps) 

head(kma_data)    # sum_rn : 일강수량 (mm)
table(kma_data$stnlds) 
length(table(kma_data$stnlds)) # 102개 지점
table(substr(kma_data$time,1,4))
length(table(substr(kma_data$time,1,4))) # 115개 년도

kma_data <- bind_cols(ID = 1:nrow(kma_data), kma_data) %>% as_tibble()
kma_data_NA_rm <- kma_data %>% filter(!is.na(sum_rn)) %>% mutate(obsyear = as.numeric(substr(time,1,4))) %>% mutate(pr = sum_rn)

### 지점별 연최대강수량을 뽑자
Pr_data <- kma_data_NA_rm %>% group_by(stnlds,obsyear,lat,long) %>% summarise(pr=max(sum_rn,na.rm=TRUE)) %>% 
  ungroup %>% mutate(obsyear = as.numeric(obsyear))


Pr_data_id <- Pr_data %>% left_join(kma_data_NA_rm, by = c('stnlds','obsyear','long','lat','pr')) %>% select(ID, stnlds, obsyear, pr, lat, long)




### 지점별로 몇년도부터 데이터가 있는지 알고싶다
table(Pr_data_id$stnlds) # 102개 지점
plot(Pr_data_id$obsyear,Pr_data_id$stnlds)

range(Pr_data_id[Pr_data_id$stnlds==285,]$obsyear)
### 1973년부터 2018년 : 46개년
stn1 <- names(table(Pr_data_id$stnlds))[c(unname(which(table(Pr_data_id$stnlds)>=46)))] %>% as.numeric() # 46개년 이상인 지점번호
Pr_46 <- Pr_data_id %>% filter(stnlds %in% stn1 , obsyear >= 1973)
length(table(Pr_46$stnlds)) # 60개 지점
# Pr_46 %>% group_by(stnlds,lat,long) %>% summarise(minyear=min(as.numeric(obsyear)), maxyear=max(as.numeric(obsyear))) %>% View


# 제주도 3개 지점, 울릉도 1개 지점 제외
unique(Pr_46[Pr_46$lat < 34,]$stnlds)
unique(Pr_46[Pr_46$long > 130,]$stnlds)
Pr_46 <- Pr_46 %>% filter(!stnlds %in% c(184,188,189,115))

# save(kma_data, file="kma_data.RData")
# save(Pr_46, file="Pr_46.RData")

