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

Do this for carbs with their upper outliers