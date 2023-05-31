## This script will allow the aggregator module to POST files to Amplicon Repository
import requests


def post_package(fp, data, user):
    '''
    Creates a package to post to Amplicon Repository

    fp --> filepath to aggregated tar gz file
    data --> a dict containing additional information:
            {'project_name': string
            'description': string
            'private': boolean
            'project_members': list}
    '''
    homepage = 'http://host.docker.internal:8000/'
    session = requests.Session()
    cookie = session.get(homepage)
    csrf_token = cookie.cookies.get_dict()['csrftoken']

    headers = {
        'Content-Type':'multipart/form-data',
        'X-CSRFToken':csrf_token,
    }

    files = {'file': open(fp, 'rb')}
    upload_url = 'http://host.docker.internal:8000/upload_api/'

    response = requests.post(upload_url,data = data, files = files)
    print(response.status_code)



if __name__ == "__main__":
    ## to do one aggregation
    # python3 /module/src/AmpliconSuiteAggregator.py -flist /module/gpunit/inputs/input_list.txt

    homepage = 'http://127.0.0.1:8000/'
    ## get cookie first
    session = requests.Session()
    cookie = session.get(homepage)
    csrf_token = cookie.cookies.get_dict()['csrftoken']
    print(csrf_token)


    upload_file = 'aggregated.tar.gz'

    upload_url = 'http://127.0.0.1:8000/upload_api/'

    headers = {
        'Content-Type':'multipart/form-data',
        'X-CSRFToken':csrf_token,
    }
    data = {'project_name': 'testPOST Private Project2',
            'description':'a test POST request to API',
            'private':True,
            'project_members':['edwin5588']}
    files = {'file': open(upload_file, 'rb')}

    response = requests.post(upload_url,data = data, files = files)
    print(response.text)

