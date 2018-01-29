## Run these commands to install and start a vnc server on a red Hat 7 instance

sudo yum groupinstall 'Server with GUI'
sudo yum install tigervnc-server

sudo cp /lib/systemd/system/vncserver@.service /etc/systemd/system/vncserver@:1.service

#Modify <USER> with ec2-user
sudo sed -i 's/<USER>/ec2-user/g' /etc/systemd/system/vncserver@:1.service
vncserver
