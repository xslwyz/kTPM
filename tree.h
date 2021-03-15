#include "head.h"

class Tree;
class TreeNode  //��ʾ���Ľ�㣬�������ֵ�
{
    friend Tree;
public:
    int data;//����
    TreeNode * left;
    TreeNode * right;  //��Ů���ֵ�ָ��
    TreeNode(int value = 0, TreeNode* l = NULL, TreeNode* r = NULL):data(value), left(l), right(r) {}//���캯��
};
class Tree
{
public:
    Tree() { root = current = NULL; }; //������
    void setRoot(int value);
    TreeNode* getRoot() ;
    int toLeft() ;
    int toRight() ;
    int Parent() ;
    int Find(int target) ;
private:
    TreeNode* root, * current;    //��ʾ���͵�ǰָ��
    //void PreOrder(ostream& out, TreeNode* p)
    int Find(TreeNode* p, int target) ;
    //void RemovesubTree(TreeNode* p)
    int FindParent(TreeNode* t, TreeNode* p) ;

};

void Tree::setRoot(int value)
{
    root = current = new TreeNode(value);  //���������
}

TreeNode* Tree::getRoot()  //�������Ҹ�,ʹ֮��Ϊ��ǰ���
{
    if (root == NULL) { current = NULL; return NULL; }
    else { current = root; return root; }
}

int Tree ::toLeft()
{ //Ѱ�ҵ�ǰ���ĵ�һ����Ů.������Ů,��������0,��ǰָ��ΪNULL;
  //���򷵻�1,��ǰָ���Ƶ���ǰ���ĵ�һ����Ů;
    if (current != NULL && current->left != NULL)
    {
        current = current->left; return 1;
    }
    current = NULL; return 0;
}


int Tree ::toRight()
{
    //Ѱ�ҵ�ǰ������һ���ֵ�.������һ���ֵ�,��������0,��ǰָ��ΪNULL;���򷵻�1,��ǰָ���ƶ�����ǰ������һ���ֵ�.
    if (current != NULL && current->right != NULL)
    {
        current = current->right; return 1;
    }
    current = NULL; return 0;
}

int Tree ::Parent()
{//Ѱ�ҵ�ǰ����˫�׽��.����˫�׽��,��������0,��ǰָ��ΪNULL;���򷵻�1,��ǰָ���Ƶ���ǰ����˫�׽��
    TreeNode* p = current;//p���浱ǰָ��
    TreeNode* t;  
    if (current == NULL || current == root) { current = NULL; return 0; }//û��˫��
    t= root;   //t�Ӹ���㿪ʼ
    int  k = FindParent(t, p);  //�Ӹ����t��������ҽ��p��˫�׽��
    return k;
}

//t����������ָ��root
int Tree ::FindParent(TreeNode* t, TreeNode* q)
{//˽�к���:�Ӹ�ָ��t��ָ���������,�ҽ��p��˫�׽��,���� current�� 
    TreeNode* p = t->left;//�ҵ�һ������
    while (q != NULL && q != p)       //q==NULL,������;q==p,�ҵ�˫��
    {
        int i = FindParent(q, p);
        if (i!= 0) return i;
        q = q->right; //����һ������
    }
    if (q != NULL && q == p) { current = t; return 1; } //�ɹ�����1
    else return 0;    //���򷵻�0
}

int Tree ::Find(int target)
{//�������������ϸ���ֵtarget�Ľ��.�����ɹ�,��������1,���򷵻�0
    if (root == NULL) return 0;      //Ϊ����
    return Find(root, target);    //����˽�к���Find����
}

int Tree ::Find(TreeNode* p, int target)
{ //˽�к���:���ɸ�ָ��p��ָ���������������������ݳ�Ա����target�Ľ��.�����ɹ�,��������1,ͬʱ�ý���Ϊ��ǰ���;����������0
    int result = 0;
    if (p->data == target) { result = 1; current = p; } //�����ɹ�,p��Ϊ��ǰ���
    else           //����,��������
    {
        TreeNode* q = p->left;
        while (q != NULL && !(result = Find(q, target))) q = q->right;
    }
    return result;
}
