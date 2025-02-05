#!/bin/bash
# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================

#echo $*

# install-kfg.sh
#   copies file $1 from $2 to $3.
#   Will create a backup copy if the file's md5 is not found in $4 (assumed to be a user's edit)
#

FILE_NAME=$1
SRC_DIR=$2
KONFIG_DIR=$3
MD5SUMS=$SRC_DIR/kfgsums

SRC_FILE=$2/$1
TGT_FILE=$3/$1

mkdir -p ${KONFIG_DIR}

function get_md5 ()
{
    local MD5="$(which md5sum)"
    if [ "$MD5" != "" ] && [ -x "$MD5" ]
    then
        "$MD5" "$1" | awk '{print $1;}'
    else
        MD5=/usr/bin/md5sum
        if [ -x "$MD5" ]
        then
            "$MD5" "$1" | awk '{print $1;}'
        else
            MD5="$(which md5)"
            if [ "$MD5" != "" ] && [ -x "$MD5" ]
            then
                "$MD5" -q "$1"
            else
                MD5=/sbin/md5
                if [ -x "$MD5" ]
                then
                    "$MD5" -q "$1"
                else
                    echo "failed to locate md5 tool" 2>&1 1>/dev/null
                    exit 5
                fi
            fi
        fi
    fi
}

#echo "installing $1 from $2 to $3, mdsums = $4"

# create a backup if installed file has been modified by user
if [ -f ${TGT_FILE} ]
then
    md5=$(get_md5 ${TGT_FILE})
    #echo "$1 md5=$md5"
    if [ "$(grep ${md5} ${MD5SUMS})" == "" ]
    then
        # not a known version of the file; create a backup copy
        mv -b -v ${TGT_FILE} ${TGT_FILE}.orig
    fi
fi

# copy to the install location
cp ${SRC_FILE} ${TGT_FILE}
