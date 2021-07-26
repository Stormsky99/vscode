git init
git remote set-url origin https://github.com/Stormsky99/vscode.git

git config --global user.name "Stromsky99"
git config --global user.email "446554759@qq.com"

git remote add origin https://github.com/Stormsky99/vscode.git

git add test.sh

git commit -m "test"

git show-ref

git branch -M main
git push -u origin master