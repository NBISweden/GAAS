#!/bin/bash

here=$(realpath "$(dirname "$0")")

for s in /mnt/*; do
	 if [ ! -d "$s" ] ; then
	     # Ignore non-directories
	     continue
	 fi

	 if [ ! -e "$s/species" ] ; then
	     # We need a species file to run here
	     continue
	 fi

	 species=$(cat "$s/species")
	 password=$(cat "$here/apolloadminpassword")
	 port=8888
	 catalinaopts=-Xmx14G

	 if [ -e "$s/port" ]; then
	     port=$(cat "$s/port")
	 fi
	 if [ -e "$s/apolloadminpassword" ]; then
	     password=$(cat "$s/apolloadminpassword")
	 fi
	 if [ -e "$s/catalinaopts" ]; then
	     catalinaopts=$(cat "$s/catalinaopts")
	 fi

	 mkdir -p "$s/data"
	 mkdir -p "$s/postgres-data"

	 docker run --restart unless-stopped --name "webapollo-$species" -d -e CATALINA_OPTS="$catalinaopts" -e APOLLO_ADMIN_EMAIL=admin@nbis.se -e APOLLO_ADMIN_PASSWORD="$password" -v "$s/data/":/data -v "$s/postgres-data/":/var/lib/postgresql -p "$port":8080 gmod/apollo:latest
done
