#### Git ####

////////////////////Wow//////////////////////


#### Basic Workflow ####
#Git is a Version control system tool
    #Keeps track of changes to files
    #Notices conflicts between changes made by different people
    #Synchronizes files between different computers

#Git stores information in the .git directory located at the root directory of the repository

#Check the status of the repository with
git status
    #Git status shows you which files are in the staging area and files that haven't been put there yet

#Git uses a staging area to store files with changes that you want to save but haven't saved yet
#The staging area is like a box and then committing those changes ships out the box

#Look at all the changes to the repository with
git diff
    #or input a directory to only see changes to that directory
    git diff directory
#in git diff display, a means 'first version' and b means 'second version'
    #--- means lines were removed and +++ means lines were added
    #@@ tells you where the changes are being made 
        #The pairs of numbers seperated by commas tells you (x,y) The start line (x) and number of lines (y) that changed

#You commit changes to a git repository in 2 steps
    #Add file(s)
    #Commit staging area
#Add a file to the staging area with
git add filename

#Compare your files with files in the staging area
git diff -r HEAD
    #-r flag means 'compare to a particular revision'
    #HEAD is a shortcut meaning 'the most recent commit'
#Compare your file to a particular file
git diff -r HEAD path/to/file

#Commit changes
git commit -m "Your message"
    #If you accidentally mistype a commit message, you can change it using the --amend flag
    git commit --amend -m "new message"
    #Make a long message in nano by omitting the -m flag. It will open up nano.
    git commit

#View logs of the repository's history (use space bar to page down and q to quit)
git log
    #View logs of a specific file or directory
    git log path




#### Repositories ####


#3-level structure for storing information
    #Metadata contained in the commit
    #Trees track names and locations in the repository when the commit happened
    #For each file in the tree, there is a blob which contains a snapshot of the file in the tree when the commit happened

#Each commit has a unique hash. You can view specific commits using (you only need first 6 characters for the hash)
git show hash

#A hash is like an absolute path but you can also use relative paths with
HEAD
HEAD~1 #second  last
HEAD~2 #third to last etc.

#View who made changes to a file
git annotate file

#View the difference betweek two commits
git diff abc123..def456 #absolute hash
git diff HEAD..HEAD~1 #relative path

#Add files
git add file

#Tell git what to ignore by creating a file in the root directory of the repository called
.gitignore
    #Use wildcards to help '*.py' or '*.txt'

#Deleting files
git clean -n
    #Shows a list of files in the repository but not tracking
git clean -f 
    #Deletes those files... for good be careful

#Configuration
git config --list --system --global --local
    #Choose one of system (all users on computer), global (per-user) or local (per-project)
#Always make sure to set your name and email
git config --global user.name JoryCurry
git config --global user.email jorycurry1@cmail.carleton.com



#### Undo ####


#If you have changed multiple files consider making seperate commits with seperate messages
git status #To see what was modified
git add path/to/file1 #to stage one file at a time
    #If you stage the wrong file use
    git reset HEAD
git commit -m "Changed this file first"
git add path/to/file2
git commit -m "Changed this file second"

#Stage the same file multiple times before committing to periodically save your work before committing

#If you want to undo changes made to an unstaged file use
git checkout -- filename     #This discards changes forever

#If you want to undo changes made to a staged file use use both lines
git reset HEAD path/to/file     #unstage the most recent version of file
git checkout -- path/to/file

#Restore an old version of a file
git log -10 file         #To show old verion of file and its hash you wanna restore to (list of 10)
git checkout hash file
git checkout -- data        #Will restore all files in the data directory to their previous state
git checkout -- .           #Revert all files in a directory to their previous version

#Undo all your changes and unstage
git reset       #Will unstage everything
git reset HEAD data         #will unsstage any files in the data directory
git reser HEAD path/to/file         #unstage a file and its changes



#### Branches ####

#Changes done on seperate branches stay seperate unit you merge them together
#By default every repository has a master branches
#List all branches in a repository
git branch

#You can use gif diff to look at differences between branches
git diff branch_1..branch_2

#Switch Branches with
git checkout branch_name

#Remove a file with 
git rm

#Create a branch
git checkout -b branch_name

#Merge Branches
git merge source destination
    #If any conflicts come up look at them with 
    git status




#### Collaborating ####


#Make a branch new repository
git init project-name
    #Project-name is the name of the repository's new root directory

#Make an existing directory a git repository
#Navigate to directory
git init
git status

#Do not create a repository inside another repository!!!!

#Clone a repository
git clone URL
git clone path/to/existing/repository new_name_MuchWow

#When you clone a repository, git remembers where the original repository was
#While in a reposity view the remotes
git remote -v

#Git automatically names the original remote 'origin' but you can add other remotes
git remote add remote-name URL
    #Remove remotes
    git remote rm remote-name
#Adding remotes allows you connect two repositories together

#Typically the workflow starts with you pulling in all up-to-date files and information from a remote
#Pulling in changes from a remote repository
git pull remote branch

#Send your changes back to the remote repository to update it with your changes
git push remote-name branch-name
#Before executing a push, always do a pull to make sure you have the most up to date stuff... and so you don't overwrite someone elses work




