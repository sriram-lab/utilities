import requests

url = "http://api.xialab.ca/mapcompounds"

payload = "{\n\t\"queryList\": \"1,3-Diaminopropane;2-Ketobutyric acid;2-Hydroxybutyric acid;\",\n\t\"inputType\": \"name\"\n}"
headers = {
    'Content-Type': "application/json",
    'cache-control': "no-cache",
    }

response = requests.request("POST", url, data=payload, headers=headers)

print(response.text)
