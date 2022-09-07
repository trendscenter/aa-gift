FROM trendscenter/gift:v1.0.0
COPY . /app

COPY ./groupicatv4.0c/icatb/nipype-0.10.0/nipype/interfaces/gift /usr/local/lib/python3.6/site-packages/nipype/interfaces/gift
