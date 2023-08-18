#!/bin/bash

rm -rf ./R/*
rm -rf ./restart_R/*
rm -rf ./restart_V/*
rm -rf ./Vk/*
rm -rf ./y_k/*
rm -rf ./xtk/*
rm -rf ./surf/*

make

./dzeq

exit 



