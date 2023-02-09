DeepAProt is a user-friendly and open-sourced deep learning model developed at Division of Agricultural Bioinformatics, ICAR-Indian Agricultural Statistics Research Institute, New Delhi for abiotic protein sequence classification. Novel activation function, namely, Gaussian Error Linear Unit with Sigmoid (SIELU) which implemented in the Deep Learning (DL)model along with other hyper parameters for classification of unknown abiotic stress protein sequences from crops of Poaceae family.  It can be accessed using http://login1.cabgrid.res.in:5500/ for the classification. Or, one can also run locally with following steps:

1. Download Pycharm using Download PyCharm: Python IDE for Professional Developers by JetBrains  and install in the local system
2. Copy the entire folder in C:\\Users\Usr\PycharmProject
3. Install the following package using the following commands:
	a. pip install sielu
	b. pip install numpy==1.19.5
	c. pip install xlsxwriter
	d. pip install pandas
	e. pip install tensorflow==2.2.0
	f. pip install keras==2.4.3
	g. pip install sklearn
	h. pip install flask
	i. pip install Bio
	j. pip install filelock

4.	 Now run the app.py file using “python app.py” command
