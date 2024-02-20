## This script will allow the aggregator module to POST files to Amplicon Repository
import requests
import sys
import os
import uuid


def post_package(fp, data, server):
    '''
    Creates a package to post to Amplicon Repository

    fp --> filepath to aggregated tar gz file
    data --> a dict containing additional information:
            {'project_name': string,
            'description': string,
            'publication_link:string,
            'private': boolean,
            'project_members': list,
            'accept_license': boolean}
    '''
    if server == 'prod':
        homepage = 'https://ampliconrepository.org'
    elif server == 'dev':
        homepage = 'https://dev.ampliconrepository.org'
    elif server == 'local-debug':
        homepage = 'http://127.0.0.1:8000'
    else:
        sys.stderr.write("Unrecognized server option: " + server + "\n")
        sys.exit(1)
    session = requests.Session()
    cookie = session.get(homepage)
    csrf_token = cookie.cookies.get_dict()['csrftoken']
    upload_file = os.path.join(fp)
    upload_url = f'{homepage}/upload_api/'
    generate_uuid = uuid.uuid4().hex
    print(generate_uuid)


    if os.path.getsize(upload_file) > 1000000000:
        ### if file is ALMOST greater than 2GB, split it up into multiple pieces
        os.system(f'mkdir -p {generate_uuid}')
        ## split into 200mb chunks for POST
        os.system(f'split -b 200m {upload_file} ./{generate_uuid}/POST_part_')
        last_part = os.listdir(f'./{generate_uuid}')[-1]
        actual_proj_name = data['project_name']
        data['project_name'] = f'MULTIPART__{generate_uuid}__{last_part}__{actual_proj_name}'
        
        for file in os.listdir(f'./{generate_uuid}'):
            fp = os.path.join(f'./{generate_uuid}', file)
            print(fp)
            files = {'file': open(fp, 'rb')}

            response = requests.post(upload_url, data = data, files = files)
            print(response.status_code)

            os.remove(fp)
        print("Large projects may take >30 minutes to be unpacked and registered after upload to the site.")
        os.rmdir(f'./{generate_uuid}')


    else:
        files = {'file': open(upload_file, 'rb')}
        response = requests.post(upload_url, data = data, files = files)
        print(response.status_code)





# if __name__ == "__main__":
#     ## to do one aggregation
#     # python3 /module/src/AmpliconSuiteAggregator.py -flist /module/gpunit/inputs/input_list.txt

#     homepage = 'http://127.0.0.1:8000'
#     ## get cookie first
#     session = requests.Session()
#     cookie = session.get(homepage)
#     csrf_token = cookie.cookies.get_dict()['csrftoken']
#     print(csrf_token)


#     upload_file = os.path.join('/Users/edwinhuang/Desktop/AA_159_testing/aggregated_results/091123_51p3MBproj.tar.gz')
#     # upload_file = '/Users/edwinhuang/Downloads/whynot.tar.gz'
#     upload_url = f'{homepage}/upload_api/'

#     headers = {
#         # 'Content-Type':'multipart/form-data',
#         'Content-Length':'100',
#         # 'X-CSRFToken':csrf_token,
#         # 'Accept-Encoding': 'gzip, deflate, br',
#     }
#     print(os.path.getsize(upload_file))
#     print(os.path.getsize('/Users/edwinhuang/Desktop/AA_159_testing/aggregated_results/091123_1p99gb_api.tar.gz'))
#     generate_uuid = uuid.uuid4().hex
#     print(generate_uuid)
#     data = {'project_name': f'TEST_NO_API_ID',
#             'description':'a test POST request to API',
#             'private':True,
#             'project_members':['edh021@ucsd.edu']}
    
#     # files = {'file': open(upload_file, 'rb')}
    
#     ## 1. create multiple tar files, for each one, post individually to website
#     ## server side: if project name is the same and contains the keyword *MULTIPART* 
#     ## then save the file to the same directory, and create a new 


#     # req = requests.Request('POST', url = upload_url, data = data, files = files)
#     # prep = session.prepare_request(req)
#     # prep.headers.update({'Content-Length':'1000'})

#     if os.path.getsize(upload_file) > 1000000000:
#         ### if file is ALMOST greater than 1 GB
#         os.system(f'mkdir -p {generate_uuid}')
#         ## split into 200mb chunks for POST
#         os.system(f'split -b 100m {upload_file} ./{generate_uuid}/POST_part_')
#         last_part = os.listdir(f'./{generate_uuid}')[-1]
#         data['project_name'] = data['project_name'] + f"__{last_part}"
        
#         for file in os.listdir(f'./{generate_uuid}'):
#             fp = os.path.join(f'./{generate_uuid}', file)
#             print(fp)
#             files = {'file': open(fp, 'rb')}

#             response = requests.post(upload_url,data = data, files = files)
#             print(response.status_code)

#             os.remove(fp)
#     else:
#         files = {'file': open(upload_file, 'rb')}
#         response = requests.post(upload_url,data = data, files = files)
#         print(response.status_code)






