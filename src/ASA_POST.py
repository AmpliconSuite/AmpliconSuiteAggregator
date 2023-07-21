## This script will allow the aggregator module to POST files to Amplicon Repository
import requests
import sys


def post_package(fp, data, server):
    '''
    Creates a package to post to Amplicon Repository

    fp --> filepath to aggregated tar gz file
    data --> a dict containing additional information:
            {'project_name': string
            'description': string
            'private': boolean
            'project_members': list}
    '''
    if server == 'prod':
        homepage = 'https://ampliconrepository.org'
    elif server == 'dev':
        homepage = 'https://dev.ampliconrepository.org'
    else:
        sys.stderr.write("Unrecognized server option: " + server + "\n")
        sys.exit(1)
    session = requests.Session()
    cookie = session.get(homepage)
    csrf_token = cookie.cookies.get_dict()['csrftoken']

    headers = {
        'Content-Type':'multipart/form-data',
        'X-CSRFToken':csrf_token,
    }

    files = {'file': open(fp, 'rb')}
    upload_url = f'{homepage}/upload_api/'

    response = requests.post(upload_url,data = data, files = files)
    print(response.status_code)
    print(response.text)



# if __name__ == "__main__":
#     ## to do one aggregation
#     # python3 /module/src/AmpliconSuiteAggregator.py -flist /module/gpunit/inputs/input_list.txt

#     homepage = 'https://ampliconrepository.org'
#     ## get cookie first
#     session = requests.Session()
#     cookie = session.get(homepage)
#     csrf_token = cookie.cookies.get_dict()['csrftoken']
#     print(csrf_token)


#     upload_file = os.path.join('/Users/edwinhuang/Downloads/test_POST_from_gpBeta.tar.gz')

#     upload_url = f'{homepage}/upload_api/'

#     headers = {
#         'Content-Type':'multipart/form-data',
#         'X-CSRFToken':csrf_token,
#     }
    
#     data = {'project_name': 'testPOST_Private_Project_with_s3adjustment_to_prod',
#             'description':'a test POST request to API',
#             'private':True,
#             'project_members':['edh021@ucsd.edu']}
    
#     files = {'file': open(upload_file, 'rb')}

#     response = requests.post(upload_url,data = data, files = files)
#     print(response.status_code)


