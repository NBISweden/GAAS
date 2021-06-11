data "external" "agent_keys" {
  program = ["bash", "-c",  "ssh-add -L | grep -v cert | jq -R -s 'split(\"\\n\")[:-1] | map({ \"key\" : . }) | .[0]'" ]
}

resource "openstack_compute_keypair_v2" "my-cloud-key" {
  name       = "${var.keypair_name}${var.project_suffix}"
  public_key =  data.external.agent_keys.result.key
}

resource "openstack_networking_network_v2" "net_webap" {
  name           = "${var.internal_network_name}${var.project_suffix}"
  admin_state_up = "true"
}

resource "openstack_networking_subnet_v2" "subnet_webap" {
  name       = "${var.internal_subnet_name}${var.project_suffix}"
  network_id = openstack_networking_network_v2.net_webap.id
  cidr       = "10.0.0.0/24"
  ip_version = 4
  dns_nameservers = ["8.8.8.8"]
  enable_dhcp = "true"
}

data "openstack_networking_network_v2" "externalnetwork" {
  name = var.external_network_name
}

resource "openstack_networking_router_v2" "router_webap" {
  name                = "${var.router_name}${var.project_suffix}"
  admin_state_up      = "true"
  external_network_id = data.openstack_networking_network_v2.externalnetwork.id
}

resource "openstack_networking_router_interface_v2" "router_interface" {
  router_id = openstack_networking_router_v2.router_webap.id
  subnet_id = openstack_networking_subnet_v2.subnet_webap.id
}



resource "openstack_compute_secgroup_v2" "secgroup_webap" {
  name        = "${var.secgroup_name}${var.project_suffix}"
  description = "security group ${var.project_suffix}"

  rule {
    from_port   = 22
    to_port     = 22
    ip_protocol = "tcp"
    cidr        = "0.0.0.0/0"
  }
  rule {
    from_port   = 443
    to_port     = 443
    ip_protocol = "tcp"
    cidr        = "0.0.0.0/0"
  }

  rule {
    from_port   = 80
    to_port     = 80
    ip_protocol = "tcp"
    cidr        = "0.0.0.0/0"
  }
  rule {
    from_port   = 8080
    to_port     = 8080
    ip_protocol = "tcp"
    cidr        = "0.0.0.0/0"
  }
  rule {
    from_port   = 8888
    to_port     = 8888
    ip_protocol = "tcp"
    cidr        = "0.0.0.0/0"
  }
}

resource "openstack_networking_port_v2" "port_webap" {
  name               = "${var.port_name}${var.project_suffix}"
  network_id         = openstack_networking_network_v2.net_webap.id
  admin_state_up     = "true"
  security_group_ids = [openstack_compute_secgroup_v2.secgroup_webap.id]

  fixed_ip {
    subnet_id  = openstack_networking_subnet_v2.subnet_webap.id
    ip_address = "10.0.0.10"
  }
}

resource "openstack_networking_floatingip_v2" "ip_webap" {
  pool = var.external_network_name
}

resource "openstack_blockstorage_volume_v2" "vol_webap" {
  name = "vol_webap${var.project_suffix}"
  size = var.size_volume
}




resource "openstack_compute_instance_v2" "instance_webap" {
  depends_on = [openstack_networking_subnet_v2.subnet_webap]
  name            = "${var.instance_name}${var.project_suffix}"
  image_name      = var.image_name
  flavor_name     = var.flavor_name
  key_pair        = openstack_compute_keypair_v2.my-cloud-key.id
  security_groups = [openstack_compute_secgroup_v2.secgroup_webap.id ]

  network {
    name = "${var.internal_network_name}${var.project_suffix}"
  }

}

resource "openstack_compute_volume_attach_v2" "volume_attachement_1" {
  instance_id = openstack_compute_instance_v2.instance_webap.id
  volume_id   = openstack_blockstorage_volume_v2.vol_webap.id
}


resource "openstack_compute_floatingip_associate_v2" "ip_webap" {
  floating_ip = openstack_networking_floatingip_v2.ip_webap.address
  instance_id = openstack_compute_instance_v2.instance_webap.id
  fixed_ip    = openstack_compute_instance_v2.instance_webap.network[0].fixed_ip_v4
}

resource "null_resource" "provision" {
  depends_on = [openstack_compute_floatingip_associate_v2.ip_webap]
  connection {
    type = "ssh"
    user = var.ssh_user
    agent  = "true"
    timeout  = "5m"
    host = openstack_networking_floatingip_v2.ip_webap.address
  }

  provisioner "local-exec" {
    command  =  "bash localscripts/get-sshkeys.sh"
  }
  provisioner "local-exec" {
    command  =  "bash localscripts/get-gaas.sh"
  }

provisioner "file" {
    source      = "tmp/annotation-cluster/ansible-ubuntu-18.04/roles/common/files/authorized-keys"
    destination = "/home/${var.ssh_user}"
  }

provisioner "remote-exec" {
  inline = [ "mkdir -p /home/${var.ssh_user}/nbis",
  	   "mkdir -p /home/${var.ssh_user}/setup",
  ]
}

provisioner "file" {
    source      = "tmp/GAAS"
    destination = "/home/${var.ssh_user}/nbis"
  }

provisioner "file" {
    source      = "remotescripts/"
    destination = "/home/${var.ssh_user}/setup"
  }

provisioner "remote-exec" {
    inline = [
      "set -x -e",
      "PATH=$PATH:/snap/bin",
      "sudo apt update",
      "sudo DEBIAN_FRONTEND=noninteractive apt dist-upgrade -y",
      "sudo DEBIAN_FRONTEND=noninteractive apt install -y ruby",
      "bash $HOME/setup/add-sshkeys.sh",
      "bash $HOME/setup/fix_repo_permissions.sh",
      "bash $HOME/setup/mount_volume.sh ${var.project_suffix}",
      "bash $HOME/setup/install_docker.sh",
      "bash $HOME/setup/install_conda.sh",
      "echo ${var.admin_password} > ~/nbis/GAAS/annotation/tools/webapollo/apollo/apolloadminpassword",
      "echo 'APOLLO_DATA_DIR=$(ls -1d /mnt/*/data | head -1)' >> $HOME/.bashrc",
      "echo export APOLLO_DATA_DIR >> $HOME/.bashrc"
    ]
}

# Do a new connection to pick up groups.

provisioner "remote-exec" {
    inline = [
          "bash $HOME/nbis/GAAS/annotation/tools/webapollo/apollo/run_webapollo.sh"
]
}

}
output "connection-string" {
    description = "To run, connect to"
    value = "${var.ssh_user}@${openstack_networking_floatingip_v2.ip_webap.address}"
}

