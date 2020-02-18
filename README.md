# Anomaly detection

Anomaly detection on a production line using principal component analysis (PCA) and kernel principal component analysis (KPCA).

## Quick preview

- _Authors_: Sylvain Combettes, Houssam L'Ghoul
- _Date_: Oct. 2018 - June 2019
- _Context_: For our penultimate-year project at Mines Nancy (half a day per week), we did research for the French company [Saint-Gobain](https://www.saint-gobain.com/en), the European or worldwide leader in all of its businesses (mainly construction materials). In 2018, Saint-Gobain had a €41.8 billion turnover, operated in 67 countries and had more than 180,000 employees.
- _Topic_: Detection of sensor failure in a production line.
- _Methods_: Principal component analysis (PCA) and kernel principal component analysis (KPCA).
- _Programming_: MATLAB.
- _Result_: the algorithm can detect 100% of the failure days observed by Saint-Gobain.
- _Links_: [report incoming]

## How to use this repository

- `datav3.mat` is a file containing data without anomalies
- `dataDefautv3.mat` is a file containing data with anomalies
- `ACP_lineaire_cstr.m` is a MATLAB script detecting anomalies in `dataDefautv3.mat` with comparison to `datav3.mat` using a (linear) PCA (principal component analysis)
- `ACP_non_lineaire_cstr.m` is a MATLAB script detecting anomalies in `dataDefautv3.mat` with comparison to `datav3.mat` using a (non-linear) KPCA (kernel principal component analysis)

`ACP_lineaire_cstr.m` can be used `ACP_non_lineaire_cstr.m` independently: there are two methods with the same goal.

## To note

- I was only able to publish one fourth of the total project, the rest being confidential.
- The MATLAB scripts `ACP_lineaire_cstr.m` and `ACP_non_lineaire_cstr.m` are commented in French. 
- The report (uploaded soon) is in French.
- The data `datav3.mat` and `dataDefautv3.mat` is from "Seongkyu Yoon and John MacGregor. Fault diagnosis with multivariate statistical models part i : using steady state fault signatures. Journal of Process Control, 11(4) :387 – 400, 2001"
