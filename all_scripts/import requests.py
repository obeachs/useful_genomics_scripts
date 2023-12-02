import requests
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
import os
path = '/Volumes/sesame/joerecovery/genomes/tair_promoters'
files = os.listdir(path)

'''opens chrome using chromedriver and will close windows when finished - no bugs so far'''
for file in files:
    filename = path + '/' + file
    driver = webdriver.Chrome(executable_path='/Volumes/sesame/movers/chromedriver_108')
    driver.get('http://plantpan.itps.ncku.edu.tw/promoter_multiple.php')

    upload = driver.find_element(By.NAME, 'Sequence_upload')
    email = driver.find_element(By.NAME, 'email')
    submit = driver.find_element(By.NAME, 'submit')

    upload.send_keys(filename)
    email.send_keys('obeachs@tcd.ie')
    submit.send_keys(Keys.RETURN)
    os.remove(filename)



