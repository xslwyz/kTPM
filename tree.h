#include "head.h"

class Tree;
class TreeNode  //表示树的结点，左孩子右兄弟
{
    friend Tree;
public:
    int data;//数据
    TreeNode * left;
    TreeNode * right;  //子女和兄弟指针
    TreeNode(int value = 0, TreeNode* l = NULL, TreeNode* r = NULL):data(value), left(l), right(r) {}//构造函数
};
class Tree
{
public:
    Tree() { root = current = NULL; }; //建立树
    void setRoot(int value);
    TreeNode* getRoot() ;
    int toLeft() ;
    int toRight() ;
    int Parent() ;
    int Find(int target) ;
private:
    TreeNode* root, * current;    //表示根和当前指针
    //void PreOrder(ostream& out, TreeNode* p)
    int Find(TreeNode* p, int target) ;
    //void RemovesubTree(TreeNode* p)
    int FindParent(TreeNode* t, TreeNode* p) ;

};

void Tree::setRoot(int value)
{
    root = current = new TreeNode(value);  //建立根结点
}

TreeNode* Tree::getRoot()  //用来查找根,使之成为当前结点
{
    if (root == NULL) { current = NULL; return NULL; }
    else { current = root; return root; }
}

int Tree ::toLeft()
{ //寻找当前结点的第一个子女.若无子女,则函数返回0,当前指针为NULL;
  //否则返回1,当前指针移到当前结点的第一个子女;
    if (current != NULL && current->left != NULL)
    {
        current = current->left; return 1;
    }
    current = NULL; return 0;
}


int Tree ::toRight()
{
    //寻找当前结点的下一个兄弟.若无下一个兄弟,则函数返回0,当前指针为NULL;否则返回1,当前指针移动到当前结点的下一个兄弟.
    if (current != NULL && current->right != NULL)
    {
        current = current->right; return 1;
    }
    current = NULL; return 0;
}

int Tree ::Parent()
{//寻找当前结点的双亲结点.若无双亲结点,则函数返回0,当前指针为NULL;否则返回1,当前指针移到当前结点的双亲结点
    TreeNode* p = current;//p保存当前指针
    TreeNode* t;  
    if (current == NULL || current == root) { current = NULL; return 0; }//没有双亲
    t= root;   //t从根结点开始
    int  k = FindParent(t, p);  //从根结点t后序遍历找结点p的双亲结点
    return k;
}

//t传进来就是指向root
int Tree ::FindParent(TreeNode* t, TreeNode* q)
{//私有函数:从根指针t所指结点后序遍历,找结点p的双亲结点,记在 current中 
    TreeNode* p = t->left;//找第一棵子树
    while (q != NULL && q != p)       //q==NULL,无子树;q==p,找到双亲
    {
        int i = FindParent(q, p);
        if (i!= 0) return i;
        q = q->right; //找下一棵子树
    }
    if (q != NULL && q == p) { current = t; return 1; } //成功返回1
    else return 0;    //否则返回0
}

int Tree ::Find(int target)
{//在树中搜索符合给定值target的结点.搜索成功,函数返回1,否则返回0
    if (root == NULL) return 0;      //为空树
    return Find(root, target);    //调用私有函数Find搜索
}

int Tree ::Find(TreeNode* p, int target)
{ //私有函数:在由根指针p所指的树或子数中搜索其数据成员等于target的结点.搜索成功,函数返回1,同时该结点成为当前结点;否则函数返回0
    int result = 0;
    if (p->data == target) { result = 1; current = p; } //搜索成功,p成为当前结点
    else           //否则,继续搜索
    {
        TreeNode* q = p->left;
        while (q != NULL && !(result = Find(q, target))) q = q->right;
    }
    return result;
}
