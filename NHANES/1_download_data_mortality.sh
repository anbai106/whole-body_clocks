mkdir -p /Users/hao/Dropbox/NHANES/linked_mortality_2019_public

cd /Users/hao/Dropbox/NHANES/linked_mortality_2019_public

for f in \
NHANES_1999_2000_MORT_2019_PUBLIC.dat \
NHANES_2001_2002_MORT_2019_PUBLIC.dat \
NHANES_2003_2004_MORT_2019_PUBLIC.dat \
NHANES_2005_2006_MORT_2019_PUBLIC.dat \
NHANES_2007_2008_MORT_2019_PUBLIC.dat \
NHANES_2009_2010_MORT_2019_PUBLIC.dat \
NHANES_2011_2012_MORT_2019_PUBLIC.dat \
NHANES_2013_2014_MORT_2019_PUBLIC.dat \
NHANES_2015_2016_MORT_2019_PUBLIC.dat \
NHANES_2017_2018_MORT_2019_PUBLIC.dat
do
  curl -k -L --fail --retry 5 --retry-delay 3 \
    "https://ftp.cdc.gov/pub/health_statistics/NCHS/datalinkage/linked_mortality/${f}" \
    -o "${f}"
done

ls -lh NHANES_*_MORT_2019_PUBLIC.dat