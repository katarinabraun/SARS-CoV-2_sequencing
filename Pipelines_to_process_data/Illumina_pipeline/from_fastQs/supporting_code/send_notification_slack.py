import os
import slack
import time

def getUserID(username, client):
	call = client.api_call("users.list")
	assert call["ok"]
	
	userlist = call["members"]
	users = filter(lambda x: x['name']==username, userlist)
	user = list(users)[0]
	userID = user['id']

	return userID

def sendSlack(user, message, counter=0):
	try:
		client = slack.WebClient(token=SLACK_BOT_TOKEN)
		
		response = client.chat_postMessage(
			channel= getUserID(user, client),
			text=message)
		assert response["ok"]
		assert response["message"]["text"] == message
	except:
		if counter < 5:
			time.wait(60)
			counter += 1
			sendSlack(user, message, counter)
		else:
			pass