Data Cleaning Read Me:

When reading in the NIR carbohydrate output from 
"Data/RawData/wsmdp_allsamples_inclval_wsmdp2021glstpls.txt"
"Data/RawData/wsmdp_allsamples_inclval_wsmdp2021glsupls.txt"
There were a number of pieces of data that were outliers. This is my determination of what to remove from the data set


First find all samples where Sucrose was < 0 
If sucrose < 0 investigate other samples of same variety and see what status of the outlier sample is 
If outlier by grubb's test, check other carb numbers. If Only the sucrose is weird set to NA
If all the carbs are weird, remove the sample entirely

Do this for all carbs with less than 0% samples

Do this for carbs with their upper outliers. Anything above Q3 + 1.5IQR check if it was an outlier



On 2021.03.19 I went through the WSMDP_2014-2015_WINY_SampleInfo.csv file and cleaned variety names
using OpenRefine, I found varieties with slightly different spellings and standardized them
for the variety W822GSe and W822Gse, I determined that their differences in capitalization were purposeful
I double checked within the 2015 field book to make sure the capitalization was correct. And then I changed
the endo mutant denotation for W822GSe varieties from se to field since they are wild type.
Also, in checking the INBHYB documentation I found that varieties denoted with sm*sm smooth and sm sm wrinkled are probably also 
W822GSe and W822Gse rewpectively. I have changed their variety names and endo sperm denotation appropriately. 