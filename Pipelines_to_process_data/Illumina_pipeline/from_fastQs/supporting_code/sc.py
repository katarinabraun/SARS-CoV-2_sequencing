import os
from send_notification_slack import sendSlack

def threadbalancer(threads,jobtype=None):
	pass

def checkexists(path):
	path = os.path.abspath(path)
	if not os.path.isdir(path):
		os.mkdir(path)

def cpu_count(num):
	if num <=1:
		return 1,1
	elif num <= 5:
		return 1,num
	else:
		results = []
		for n in range(2,num):
			if num % n == 0:
				results.append(n)
		if len(results) < 2:
			return cpu_count(num-1)
		index = int(len(results)/2)
		return results[index],int(num/results[index])


def procTitle(title, cfg):
	printtitle = ' '+title.strip()+' '
	rows,cols = os.popen('stty size','r').read().split()
	print("\n" + printtitle.center(59,'*') + "\n")
	if cfg['slackUser']:
		sendSlack(cfg['slackUser'], title)
	with open(cfg['exec']['logfile'],'a') as outlog:
		outlog.write('***********\n')
		outlog.write(f'{title}\n')

def mainTitle():
	print('''
 #####
#     #  #    #  #  ######  ######  #       ######   ####
#        ##   #  #  #       #       #       #       #
 #####   # #  #  #  #####   #####   #       #####    ####
      #  #  # #  #  #       #       #       #            #
#     #  #   ##  #  #       #       #       #       #    #
 #####   #    #  #  #       #       ######  ######   ####
	''')
	print('\n')
