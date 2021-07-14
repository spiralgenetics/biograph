#!/usr/bin/python3
#
# Generates an AMI that can be used to start up AWS instances.
# Currently only developer AMIs are supported.
#
# See https://spiralgenetics.atlassian.net/wiki/display/DEV/Developer+AMIs for usage details.

import argparse
import urllib2
import string
import boto3.ec2
import uuid
import time
import os

image_types = {
    "dev": ["common", "dev"],
    "jenkins": ["common", "dev", "jenkins"],
}

# distribution should be one of xenial, trusty, etc.
def get_ubuntu_ami_id(distribution, region):
    ami_id = None
    release_date = None
    url = "http://cloud-images.ubuntu.com/query/" + distribution + "/server/released.current.txt"
    print "Querying: " + url
    releases = urllib2.urlopen(url)
    for line in releases.readlines():
        line = line.rstrip("\r\n")
        parts = string.split(line, "\t")
        if len(parts) != 11:
            print repr(parts)
            raise Exception("Expecting 9 fields in release file: " + line)
        # print repr(parts)
        if parts[2] != "release":
            continue
        if parts[4] != "ebs-ssd":
            continue
        if parts[6] != region:
            continue
        if parts[10] != "hvm":
            continue
        print "Matched line: " + line
        if not ami_id:
            ami_id = parts[7]
            release_date = parts[3]
            print "Using ami id " + ami_id + " released at " + release_date
    if not ami_id:
        raise Exception("Unable to find a matching AMI ID")
    return (ami_id, release_date)

def main():

    parser = argparse.ArgumentParser(description='Ubuntu developer AMI creation tool')
    parser.add_argument('-d', '--distribution', metavar='dist', default='xenial', help='Distribution name (default: xenial)')
    parser.add_argument('-r', '--region', default='us-west-2', metavar='region', help='AWS region (default: us-west-2)')
    parser.add_argument('-p', '--purpose', metavar='purpose', default='dev', help='Purpose (dev)')
    parser.add_argument('-k', '--keyname', default='deployment', help='SSH key name (default: deployment)')
    parser.add_argument('-s', '--sg', default='inbound-ssh-only', help='Security group name (default: ssh_only)')
    args = parser.parse_args()

    os.putenv('AWS_DEFAULT_REGION', args.region)

    (ubuntu_ami_id, release_date) = get_ubuntu_ami_id(args.distribution, args.region)

    #ubuntu_ami_id = "ami-7c803d1c"
    #releaes_date = "20170101"
    print "Ami ID: " + ubuntu_ami_id
    ec2 = boto3.resource('ec2', region_name=args.region)

    script_data = "#!/bin/bash\n"

    for setup_script in image_types[args.purpose]:
        with open(os.path.dirname(os.path.abspath(__file__)) + "/" + setup_script + "-image-setup.sh", "r") as common_setup:
            script_data = script_data + common_setup.read()
            script_data = script_data + "\n"

    # Shut down after image is prepared
    script_data = script_data + "\nhalt\n"

    instances = ec2.create_instances(ImageId = ubuntu_ami_id, MinCount = 1, MaxCount = 1,
                                     KeyName = args.keyname,
                                     UserData = script_data,
                                     InstanceType = 'm3.medium',
                                     SecurityGroups = [args.sg],
                                     IamInstanceProfile={ 'Name': 'Developer' }
                                 )

    instance_ids = []
    instance = None
    for i in instances:
        instance_ids.append(i)
        instance = i
    if len(instance_ids) != 1:
        raise Exception("Attempting to start 1 instance; got: " + repr(instances))

    print "Instance id: " + instance.id
    instance.create_tags(Tags = [{"Key": "Name",
                                  "Value": "generate-" + args.purpose + "-ami-" + str(uuid.uuid4())}])

    while True:
        time.sleep(10)
        instance.reload()
        print "Instance type now: " + repr(instance.state)
        if instance.state['Name'] == 'stopped':
            break

    print "Creating new AMI"
    # Create an ami from the instance
    image_name = args.purpose + "-base-" + args.distribution + "-" + time.strftime("%Y%m%d.%H%M%S")
    new_image = instance.create_image(Name=image_name,
                                      Description="Generated from " + ubuntu_ami_id)
    new_image.create_tags(Tags = [{"Key": "Name",
                                   "Value": image_name}])
    print "New ami id: " + new_image.id

    while True:
        time.sleep(10)
        new_image.reload()
        print "New image state: " + new_image.state
        if new_image.state == "available":
            break

    print "Terminating instance"
    instance.terminate()

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        pass
